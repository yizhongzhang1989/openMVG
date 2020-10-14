// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_HPP
#define OPENMVG_SFM_SFM_DATA_HPP

#include <string>
#include <fstream>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"

namespace openMVG {
namespace sfm {

class CsvReader
{
public:
    CsvReader() = delete;
    CsvReader( const std::string& filename, const char split = ',' )
    {
        split_= split;
        csvInput_.open( filename );
        if( !csvInput_ ) throw std::runtime_error( "Error csv file dict: " + filename );

        std::string header;
        getline( csvInput_, header );

        data_ = std::vector<double>(7, 0);
        last_data_ = std::vector<double>(7, 0);
    }
    ~CsvReader()
    {
        if(csvInput_) csvInput_.close();
    }

    bool readline()
    {
        if( !csvInput_ ) throw std::runtime_error(" Not Set File to Read ");
        if( csvInput_.eof() ) return false;

        std::string line;
        getline(csvInput_, line);

        std::istringstream readStr( line );
        std::string part;

        for( int i= 0;i<7;++i  )
        {
            getline( readStr, part, split_ );
            data_[i] = std::strtod( part.c_str(), nullptr );
        }
        if( last_data_[0] != 0 )
        {
            if( (data_[0] - last_data_[0]) != 5 )
                return false;
        }

        last_data_ = data_;

        return true;
    }

    std::vector<double> data_;
private:
    std::vector<double> last_data_;
    std::ifstream csvInput_;
    char split_;
};

class IMU_Dataset
{
public:
    explicit IMU_Dataset( const std::string&IMU_file_path )
    {
        CsvReader reader( IMU_file_path );
        while(reader.readline())
        {
            vec_times.push_back(reader.data_[0]);
            vec_acc.emplace_back(reader.data_[1],reader.data_[2],reader.data_[3]);
            vec_gyr.emplace_back(reader.data_[4],reader.data_[5],reader.data_[6]);
        }
    }

    void corect_time(const IndexT _time)
    {
        double dt = vec_times.back() - _time;
        corect_dt(dt);
    }

    void corect_dt( const IndexT _dt )
    {
        IndexT min_new = std::numeric_limits<IndexT>::max();
        for( auto&time:vec_times )
        {
            if( min_new == std::numeric_limits<IndexT>::max() && time > _dt )
            {
                min_new = time - _dt;
            }
            time -= _dt;
        }
        assert( min_new != std::numeric_limits<IndexT>::max() );
        update_measure(min_new);
    }

    void update_measure(const IndexT min_new)
    {
        std::vector<Eigen::Vector3d> vec_acc_new;
        std::vector<Eigen::Vector3d> vec_gyr_new;
        std::vector<IndexT> vec_times_new;

        size_t index =0;
        for( ; index < vec_times.size(); ++index )
        {
            if( vec_times[index] == min_new ) break;
        }
        assert( index < vec_times.size() );
        for( ; index < vec_times.size(); ++index )
        {
            vec_acc_new.push_back( vec_acc[index] );
            vec_gyr_new.push_back( vec_gyr[index] );
            vec_times_new.push_back( vec_times[index] );
        }
        vec_times = vec_times_new;
        vec_acc = vec_acc_new;
        vec_gyr = vec_gyr_new;
    }

    std::tuple< std::vector<IndexT>, std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d> > GetMeasure(const IndexT _t0, const IndexT _t1)
    {
        std::vector<Eigen::Vector3d> vec_acc_part;
        std::vector<Eigen::Vector3d> vec_gyr_part;
        std::vector<IndexT> vec_times_part;
        if( _t0 == _t1 )
            return std::make_tuple( vec_times_part, vec_acc_part, vec_gyr_part );
        size_t index =0;
        for(  ;index < vec_times.size(); ++index )
        {
            if( vec_times[index] > _t0 ) break;
        }

        assert( index < vec_times.size() );
        if( vec_times[index] > _t0 && index > 0 )
        {
            double k1, k2;
            k1 = static_cast<double>( _t0 - vec_times[index-1] )/static_cast<double>( vec_times[index] - vec_times[index-1] );
            k2 = static_cast<double>( vec_times[index] - _t0 )/static_cast<double>( vec_times[index] - vec_times[index-1] );
            Eigen::Vector3d acc0 = vec_acc[index-1];
            Eigen::Vector3d acc1 = vec_acc[index];
            Eigen::Vector3d gyr0 = vec_gyr[index-1];
            Eigen::Vector3d gyr1 = vec_gyr[index];

            Eigen::Vector3d acc_cur = k1 * acc0 + k2 * acc1;
            Eigen::Vector3d gyr_cur = k1 * gyr0 + k2 * gyr1;

            vec_acc_part.push_back( acc_cur );
            vec_gyr_part.push_back( gyr_cur );
            vec_times_part.push_back( _t0 );
        }

        for(  ;index < vec_times.size(); ++index )
        {
            if( vec_times[index] > _t1 ) break;
            vec_acc_part.push_back( vec_acc[index] );
            vec_gyr_part.push_back( vec_gyr[index] );
            vec_times_part.push_back( vec_times[index] );
        }
        if( vec_times[index-1] < _t1 && index < vec_times.size() )
        {

            double k1, k2;
            k1 = static_cast<double>( _t1 - vec_times[index-1] )/static_cast<double>( vec_times[index] - vec_times[index-1] );
            k2 = static_cast<double>( vec_times[index] - _t1 )/static_cast<double>( vec_times[index] - vec_times[index-1] );
            Eigen::Vector3d acc0 = vec_acc[index-1];
            Eigen::Vector3d acc1 = vec_acc[index];
            Eigen::Vector3d gyr0 = vec_gyr[index-1];
            Eigen::Vector3d gyr1 = vec_gyr[index];

            Eigen::Vector3d acc_cur = k1 * acc0 + k2 * acc1;
            Eigen::Vector3d gyr_cur = k1 * gyr0 + k2 * gyr1;

            vec_acc_part.push_back( acc_cur );
            vec_gyr_part.push_back( gyr_cur );
            vec_times_part.push_back( _t1 );
        }

        return std::make_tuple( vec_times_part, vec_acc_part, vec_gyr_part );
    }

    // TODO xinli change to map
private:
    std::vector<Eigen::Vector3d> vec_acc;
    std::vector<Eigen::Vector3d> vec_gyr;
    std::vector<IndexT> vec_times;
};

#define ACC_N 0.01
#define GYR_N 0.005
#define ACC_W 0.0002
#define GYR_W 4e-6

class IMU_InteBase
{
public:
    IMU_InteBase() = delete;
    ~IMU_InteBase() = default;

    IMU_InteBase(const IndexT _t0, const IndexT _t1)
            :
            jacobian{Eigen::Matrix<double, 15, 15>::Identity()}, covariance{Eigen::Matrix<double, 15, 15>::Zero()},
            linearized_ba_(Eigen::Vector3d(0.,0.,0.)), linearized_bg_{Eigen::Vector3d(0.,0.,0.)}, sum_dt_( _t1 - _t0 ), t0_(_t0), t1_(_t1),
            delta_p_{Eigen::Vector3d::Zero()}, delta_q_{Eigen::Quaterniond::Identity()}, delta_v_{Eigen::Vector3d::Zero()}

    {
        {
            sum_dt_ = (static_cast<double>(t1_) - static_cast<double>(t0_)) / 1000.;
        }
        noise = Eigen::Matrix<double, 18, 18>::Zero();
        noise.block<3, 3>(0, 0) =  (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
        noise.block<3, 3>(3, 3) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
        noise.block<3, 3>(6, 6) =  (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
        noise.block<3, 3>(9, 9) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
        noise.block<3, 3>(12, 12) =  (ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
        noise.block<3, 3>(15, 15) =  (GYR_W * GYR_W) * Eigen::Matrix3d::Identity();
    }

    void integrate( const std::vector< Eigen::Vector3d >& _accs, const std::vector<Eigen::Vector3d>& _gyrs, const std::vector<IndexT>& _times_T )
    {
        double last_t = static_cast<double>(t0_);
        last_t /= 1000.;
        linearized_acc_ = _accs[0];
        linearized_gyr_ = _gyrs[0];
        for( size_t index = 0; index < _times_T.size(); ++index )
        {
            double time = static_cast<double>(_times_T[index]);
            time /= 1000.;
            double dt = time - last_t;
            last_t = time;
            propagate(dt, _accs[index], _gyrs[index]);
            dt_buf_.push_back(dt);
        }
        acc_buf_ = _accs;
        gyr_buf_ = _gyrs;
    }

    void repropagate(const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg)
    {
        sum_dt_ = 0.0;
        acc_0_ = linearized_acc_;
        gyr_0_ = linearized_gyr_;
        delta_p_.setZero();
        delta_v_.setZero();
        delta_q_.setIdentity();
        linearized_ba_ = _linearized_ba;
        linearized_bg_ = _linearized_bg;
        jacobian.setIdentity();
        covariance.setZero();
        for (int i = 0; i < static_cast<int>(dt_buf_.size()); i++)
            propagate(dt_buf_[i], acc_buf_[i], gyr_buf_[i]);
    }

    void propagate(const double dt, const Eigen::Vector3d& acc_1, const Eigen::Vector3d& gyro_1)
    {
        acc_1_ = acc_1;
        gyr_1_ = gyro_1;

        Eigen::Vector3d result_delta_p;
        Eigen::Quaterniond result_delta_q;
        Eigen::Vector3d result_delta_v;
        Eigen::Vector3d result_linearized_ba;
        Eigen::Vector3d result_linearized_bg;

        if( dt > 0 )
            midPointIntegration(dt, acc_0_, gyr_0_, acc_1, gyro_1, delta_p_, delta_q_, delta_v_,linearized_ba_, linearized_bg_,
                            result_delta_p, result_delta_q, result_delta_v, result_linearized_ba, result_linearized_bg, true);

        //checkJacobian(_dt, acc_0, gyr_0, acc_1, gyr_1, delta_p, delta_q, delta_v,
        //                    linearized_ba, linearized_bg);
        delta_p_ = result_delta_p;
        delta_q_ = result_delta_q;
        delta_v_ = result_delta_v;
        sum_dt_ += dt;
        acc_0_ = acc_1;
        gyr_0_ = gyro_1;
    }

    void midPointIntegration(const double dt,
                             const Eigen::Vector3d& acc_0, const Eigen::Vector3d& gyro_0,
                             const Eigen::Vector3d& acc_1, const Eigen::Vector3d& gyro_1,
                             const Eigen::Vector3d& delta_p, const Eigen::Quaterniond& delta_q, const Eigen::Vector3d& delta_v,
                             const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg,
                             Eigen::Vector3d& result_delta_p, Eigen::Quaterniond& result_delta_q, Eigen::Vector3d& result_delta_v,
                             Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg,
                             const bool update_jacobian = false
    )
    {
        Eigen::Vector3d un_acc0 = delta_q * ( acc_0 - linearized_ba );
        Eigen::Vector3d un_gyr = .5 * (gyro_0 + gyro_1) - linearized_bg;
        result_delta_q
                = delta_q * Eigen::Quaterniond(1, un_gyr(0) * dt / 2, un_gyr(1) * dt / 2, un_gyr(2) * dt / 2);
        result_delta_q.normalize();
        Eigen::Vector3d un_acc_1 = result_delta_q * (acc_1 - linearized_ba);
        Eigen::Vector3d un_acc = .5 * (un_acc0 + un_acc_1);
        result_delta_p = delta_p + delta_v * dt + 0.5 * un_acc * dt * dt;
        result_delta_v = delta_v + un_acc * dt;
        result_linearized_ba = linearized_ba;
        result_linearized_bg = linearized_bg;


        if(update_jacobian)
        {
            Eigen::Vector3d w_x = 0.5 * (gyro_0 + gyro_1) - linearized_bg;
            Eigen::Vector3d a_0_x = acc_0 - linearized_ba;
            Eigen::Vector3d a_1_x = acc_1 - linearized_ba;
            Eigen::Matrix3d R_w_x, R_a_0_x, R_a_1_x;

            R_w_x<<0, -w_x(2), w_x(1),
                    w_x(2), 0, -w_x(0),
                    -w_x(1), w_x(0), 0;
            R_a_0_x<<0, -a_0_x(2), a_0_x(1),
                    a_0_x(2), 0, -a_0_x(0),
                    -a_0_x(1), a_0_x(0), 0;
            R_a_1_x<<0, -a_1_x(2), a_1_x(1),
                    a_1_x(2), 0, -a_1_x(0),
                    -a_1_x(1), a_1_x(0), 0;

            Eigen::MatrixXd F = Eigen::MatrixXd::Zero(15, 15);
            F.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
            F.block<3, 3>(0, 3) = -0.25 * delta_q.toRotationMatrix() * R_a_0_x * dt * dt +
                                  -0.25 * result_delta_q.toRotationMatrix() * R_a_1_x * (Eigen::Matrix3d::Identity() - R_w_x * dt) * dt * dt;
            F.block<3, 3>(0, 6) = Eigen::MatrixXd::Identity(3,3) * dt;
            F.block<3, 3>(0, 9) = -0.25 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * dt * dt;
            F.block<3, 3>(0, 12) = -0.25 * result_delta_q.toRotationMatrix() * R_a_1_x * dt * dt * -dt;
            F.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() - R_w_x * dt;
            F.block<3, 3>(3, 12) = -1.0 * Eigen::MatrixXd::Identity(3,3) * dt;
            F.block<3, 3>(6, 3) = -0.5 * delta_q.toRotationMatrix() * R_a_0_x * dt +
                                  -0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * (Eigen::Matrix3d::Identity() - R_w_x * dt) * dt;
            F.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity();
            F.block<3, 3>(6, 9) = -0.5 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * dt;
            F.block<3, 3>(6, 12) = -0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * dt * -dt;
            F.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity();
            F.block<3, 3>(12, 12) = Eigen::Matrix3d::Identity();
            //cout<<"A"<<endl<<A<<endl;

            Eigen::MatrixXd V = Eigen::MatrixXd::Zero(15,18);
            V.block<3, 3>(0, 0) =  0.25 * delta_q.toRotationMatrix() * dt * dt;
            V.block<3, 3>(0, 3) =  0.25 * -result_delta_q.toRotationMatrix() * R_a_1_x  * dt * dt * 0.5 * dt;
            V.block<3, 3>(0, 6) =  0.25 * result_delta_q.toRotationMatrix() * dt * dt;
            V.block<3, 3>(0, 9) =  V.block<3, 3>(0, 3);
            V.block<3, 3>(3, 3) =  0.5 * Eigen::MatrixXd::Identity(3,3) * dt;
            V.block<3, 3>(3, 9) =  0.5 * Eigen::MatrixXd::Identity(3,3) * dt;
            V.block<3, 3>(6, 0) =  0.5 * delta_q.toRotationMatrix() * dt;
            V.block<3, 3>(6, 3) =  0.5 * -result_delta_q.toRotationMatrix() * R_a_1_x  * dt * 0.5 * dt;
            V.block<3, 3>(6, 6) =  0.5 * result_delta_q.toRotationMatrix() * dt;
            V.block<3, 3>(6, 9) =  V.block<3, 3>(6, 3);
            V.block<3, 3>(9, 12) = Eigen::MatrixXd::Identity(3,3) * dt;
            V.block<3, 3>(12, 15) = Eigen::MatrixXd::Identity(3,3) * dt;

            //step_jacobian = F;
            //step_V = V;
            jacobian = F * jacobian;
            covariance = F * covariance * F.transpose() + V * noise * V.transpose();
        }
    }

    double sum_dt_; // scond
    IndexT t0_;  // second * 1000
    IndexT t1_;  // second * 1000

private:

    std::vector<double> dt_buf_;
    std::vector<Eigen::Vector3d> acc_buf_;
    std::vector<Eigen::Vector3d> gyr_buf_;
    Eigen::Vector3d linearized_ba_, linearized_bg_;
    Eigen::Vector3d linearized_acc_, linearized_gyr_;

    Eigen::Vector3d acc_0_, gyr_0_;
    Eigen::Vector3d acc_1_, gyr_1_;

    Eigen::Vector3d delta_p_;
    Eigen::Vector3d delta_v_;
    Eigen::Quaterniond delta_q_;

    Eigen::Matrix<double, 15, 15> jacobian, covariance;
    Eigen::Matrix<double, 18, 18> noise;
};

using Imus = Hash_Map<IndexT, IMU_InteBase>;

using Timestamps = Hash_Map<IndexT, double>;

/// Define a collection of IntrinsicParameter (indexed by View::id_intrinsic)
using Intrinsics = Hash_Map<IndexT, std::shared_ptr<cameras::IntrinsicBase>>;

/// Define a collection of Pose (indexed by View::id_pose)
using Poses = Hash_Map<IndexT, geometry::Pose3>;

/// Define a collection of View (indexed by View::id_view)
using Views = Hash_Map<IndexT, std::shared_ptr<View>>;

/// Generic SfM data container
/// Store structure and camera properties:
struct SfM_Data
{
  /// Considered views
  Views views;
  /// Considered poses (indexed by view.id_pose)
  Poses poses;

  Imus imus;

  // indexed by view.id_pose
  Timestamps timestamps;

  std::shared_ptr<IMU_Dataset> imu_dataset;
  /// Considered camera intrinsics (indexed by view.id_intrinsic)
  Intrinsics intrinsics;
  /// Structure (3D points with their 2D observations)
  Landmarks structure;
  /// Controls points (stored as Landmarks (id_feat has no meaning here))
  Landmarks control_points;

  /// Root Views path
  std::string s_root_path;

  //--
  // Accessors
  //--
  const Views & GetViews() const {return views;}
  const Poses & GetPoses() const {return poses;}
  const Intrinsics & GetIntrinsics() const {return intrinsics;}
  const Landmarks & GetLandmarks() const {return structure;}
  const Landmarks & GetControl_Points() const {return control_points;}

  /// Check if the View have defined intrinsic and pose
  bool IsPoseAndIntrinsicDefined(const View * view) const
  {
    if (!view) return false;
    return (
      view->id_intrinsic != UndefinedIndexT &&
      view->id_pose != UndefinedIndexT &&
      intrinsics.find(view->id_intrinsic) != intrinsics.end() &&
      poses.find(view->id_pose) != poses.end());
  }

  /// Get the pose associated to a view
  const geometry::Pose3 GetPoseOrDie(const View * view) const
  {
    return poses.at(view->id_pose);
  }
};

struct IMU_Data
{
    Imus imus;
    Timestamps timestamps;
};

struct VISfM_Data : public SfM_Data
{
    Imus imus;
    Timestamps timestamps;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_HPP
