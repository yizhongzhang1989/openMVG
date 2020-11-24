//
// Created by xin on 2020/10/18.
//

#ifndef OPENMVG_IMU_INTEBASE_HPP
#define OPENMVG_IMU_INTEBASE_HPP

#include <fstream>
#include <Eigen/Core>
//#include <eigen3/Eigen/Core>
#include <vector>
#include <Eigen/Dense>
#include <openMVG/numeric/eigen_alias_definition.hpp>
#include "Utility.hpp"
#include "VI_static_Parm.hpp"
#include "openMVG/types.hpp"
#include <iostream>

namespace openMVG
{
    namespace sfm
    {
        enum class IMU_File_Type : int
        {
            Mate20Pro                = 1,
            EuRoc     = 2
        };

        class CsvReader
        {
        public:
            CsvReader() = delete;
            CsvReader( const std::string& filename, IMU_File_Type imuFileType, const char split = ',' )
            {
                split_= split;
                imuFileType_ = imuFileType;
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

                if( imuFileType_ == IMU_File_Type::EuRoc )
                {
                    data_[0] -= 1403715 * 1e12;
                    data_[0] /= 1e6;
                    data_[0] = static_cast<long long int>(data_[0]);
                }

                if( last_data_[0] != 0 )
                {
                    if( (data_[0] - last_data_[0]) != 5 )
                    {
                        std::cout << line << std::endl;
                        std::cout << data_[0] << " - " <<last_data_[0] << " != 5" << std::endl;
                        return false;
                    }
                }

//                if( last_data_[0] == 46274 )
//                {
//                    return false;
//                }


                last_data_ = data_;

                return true;
            }

            std::vector<double> data_;
        private:
            std::vector<double> last_data_;
            std::ifstream csvInput_;
            IMU_File_Type imuFileType_;
            char split_;
        };

        class IMU_Dataset
        {
        public:
            explicit IMU_Dataset( const std::string&IMU_file_path, const std::string& simu_file_type )
            {
                IMU_File_Type imu_file_type;
                if( simu_file_type == std::string( "Mate20Pro" ) || simu_file_type == std::string( "Simu" ) )
                {
                    imu_file_type = IMU_File_Type::Mate20Pro;
                    std::cout << "imu_file_type = IMU_File_Type::Mate20Pro;" << std::endl;
                }
                else if( simu_file_type == std::string("EuRoc") )
                {
                    imu_file_type = IMU_File_Type::EuRoc;
                    std::cout << "imu_file_type = IMU_File_Type::EuRoc;" << std::endl;
                }


                CsvReader reader( IMU_file_path, imu_file_type );
                while(reader.readline())
                {
                    switch (imu_file_type) {
                        case IMU_File_Type::Mate20Pro:
                        {
                            vec_times.push_back(reader.data_[0]);
                            vec_acc.emplace_back(reader.data_[1],reader.data_[2],reader.data_[3]);
                            vec_gyr.emplace_back(reader.data_[4],reader.data_[5],reader.data_[6]);
                            break;
                        }
                        case IMU_File_Type::EuRoc:
                        {
                            vec_times.push_back(reader.data_[0]);
                            vec_gyr.emplace_back(reader.data_[1],reader.data_[2],reader.data_[3]);
                            vec_acc.emplace_back(reader.data_[4],reader.data_[5],reader.data_[6]);
                            break;
                        }
                        default:
                            break;

                    }
                }
            }

            void corect_time(const double _time)
            {
                std::cout << "vec_times.back() = "  << vec_times.back() << std::endl;
                std::cout << "_time = " << _time << std::endl;
                double dt = vec_times.back() - _time;
                std::cout << "dt = " << dt << std::endl;
                corect_dt(dt);
            }

            void corect_dt( const double _dt )
            {
                std::cout << "_dt = " << _dt << std::endl;
//                int64_t min_new = std::numeric_limits<int64_t>::max();

                for( auto&time:vec_times )
                {
//                    if( min_new == std::numeric_limits<IndexT>::max() && time > _dt )
//                    {
//                        min_new = time - _dt;
//                    }
                    time += _dt;
                }
//                assert( min_new != std::numeric_limits<IndexT>::max() );
//                update_measure(min_new);
            }

            void update_measure(const int64_t min_new)
            {
                std::vector<Vec3> vec_acc_new;
                std::vector<Vec3> vec_gyr_new;
                std::vector<double> vec_times_new;

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

            std::tuple< bool, std::vector<double>, std::vector<Vec3>, std::vector<Vec3> > GetMeasure(const double _t0, const double _t1)
            {
                bool ret_flag = true;
                std::vector<Vec3> vec_acc_part;
                std::vector<Vec3> vec_gyr_part;
                std::vector<double> vec_times_part;
                if( _t0 == _t1 )
                    return std::make_tuple( false, vec_times_part, vec_acc_part, vec_gyr_part );
                if( vec_times[0] > _t0 )
                    return std::make_tuple( false, vec_times_part, vec_acc_part, vec_gyr_part );
                size_t index =0;
                for(  ;index < vec_times.size(); ++index )
                {
                    if( vec_times[index] >= 0 ) break;
                }
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
                    Vec3 acc0 = vec_acc[index-1];
                    Vec3 acc1 = vec_acc[index];
                    Vec3 gyr0 = vec_gyr[index-1];
                    Vec3 gyr1 = vec_gyr[index];

                    Vec3 acc_cur = k1 * acc0 + k2 * acc1;
                    Vec3 gyr_cur = k1 * gyr0 + k2 * gyr1;

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
                    Vec3 acc0 = vec_acc[index-1];
                    Vec3 acc1 = vec_acc[index];
                    Vec3 gyr0 = vec_gyr[index-1];
                    Vec3 gyr1 = vec_gyr[index];

                    Vec3 acc_cur = k1 * acc0 + k2 * acc1;
                    Vec3 gyr_cur = k1 * gyr0 + k2 * gyr1;

                    vec_acc_part.push_back( acc_cur );
                    vec_gyr_part.push_back( gyr_cur );
                    vec_times_part.push_back( _t1 );
                }

                if( index == vec_times.size() && vec_times[index-1] < _t1 )
                    return std::make_tuple( false, vec_times_part, vec_acc_part, vec_gyr_part );
                else
                    return std::make_tuple( true, vec_times_part, vec_acc_part, vec_gyr_part );
            }

            // TODO xinli change to map
//        private:
            std::vector<Vec3> vec_acc;
            std::vector<Vec3> vec_gyr;
            std::vector<double> vec_times;
        };

//#define GYR_N 9.7933408260869451e-04
//#define GYR_W 1.7270393511376410e-05
//
//#define ACC_W 9.7230832140122278e-04
//#define ACC_N 1.3061437477214574e-02


//#define ACC_N 0.08
//#define GYR_N 0.004
//#define ACC_W 0.00004
//#define GYR_W 2.0e-6

        enum StateOrder
        {
            O_P = 0,
            O_R = 3,
            O_V = 6,
            O_BA = 9,
            O_BG = 12
        };

        class IMU_InteBase
        {
        public:
//            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            IMU_InteBase()
            {
                std::runtime_error("should not be here");
            }
            ~IMU_InteBase() = default;

            IMU_InteBase(const double _t0, const double _t1)
                    :
                    jacobian{Eigen::Matrix<double, 15, 15>::Identity()}, covariance{Eigen::Matrix<double, 15, 15>::Zero()},
                    linearized_ba_(Eigen::Vector3d(0.,0.,0.)), linearized_bg_{Eigen::Vector3d(0.,0.,0.)}, t0_(_t0), t1_(_t1),
                    delta_p_{Eigen::Vector3d::Zero()}, delta_q_{Eigen::Quaterniond::Identity()}, delta_v_{Eigen::Vector3d::Zero()}

            {
                {
                    sum_dt_ = (static_cast<double>(t1_) - static_cast<double>(t0_)) / 1000.;
                }
                noise = Eigen::Matrix<double, 18, 18>::Zero();
                noise.block<3, 3>(0, 0) =  (VIstaticParm::acc_n * VIstaticParm::acc_n) * Eigen::Matrix3d::Identity();
                noise.block<3, 3>(3, 3) =  (VIstaticParm::gyr_n * VIstaticParm::gyr_n) * Eigen::Matrix3d::Identity();
                noise.block<3, 3>(6, 6) =  (VIstaticParm::acc_n * VIstaticParm::acc_n) * Eigen::Matrix3d::Identity();
                noise.block<3, 3>(9, 9) =  (VIstaticParm::gyr_n * VIstaticParm::gyr_n) * Eigen::Matrix3d::Identity();
                noise.block<3, 3>(12, 12) =  (VIstaticParm::acc_w * VIstaticParm::acc_w) * Eigen::Matrix3d::Identity();
                noise.block<3, 3>(15, 15) =  (VIstaticParm::gyr_w * VIstaticParm::gyr_w) * Eigen::Matrix3d::Identity();
            }

            void integrate( const std::vector< Vec3 >& _accs, const std::vector<Vec3>& _gyrs, const std::vector<double>& _times_T, bool good_falg = true )
            {
//                std::cout << "_accs.size() = " << _accs.size() << std::endl;
//                std::cout << "good_falg = " << good_falg
//                << std::endl;
                good_to_opti_ = good_falg;
                dt_buf_.clear();
                acc_buf_.clear();
                gyr_buf_.clear();
                delta_p_.setZero();
                delta_v_.setZero();
                delta_q_.setIdentity();
                jacobian.setIdentity();
                covariance.setZero();
                if( VIstaticParm::acc_n == 0 && VIstaticParm::gyr_n == 0 && VIstaticParm::acc_w == 0  && VIstaticParm::gyr_w == 0  )
                    covariance.setIdentity();
                double last_t = static_cast<double>(t0_);
                last_t /= 1000.;
                if(!good_falg) return;
                linearized_acc_ = _accs[0];
                linearized_gyr_ = _gyrs[0];
                acc_0_ = _accs[0];
                gyr_0_ = _gyrs[0];
                for( size_t index = 1; index < _times_T.size(); ++index )
                {
                    double time = static_cast<double>(_times_T[index]);
                    time /= 1000.;
                    double dt = time - last_t;
                    last_t = time;
                    propagate(dt, _accs[index], _gyrs[index]);
                    dt_buf_.push_back(dt);
                    acc_buf_.push_back(_accs[index]);
                    gyr_buf_.push_back(_accs[index]);
                }
            }

            void repropagate(const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg)
            {
                acc_0_ = linearized_acc_;
                gyr_0_ = linearized_gyr_;
                delta_p_.setZero();
                delta_v_.setZero();
                delta_q_.setIdentity();
                linearized_ba_ = _linearized_ba;
                linearized_bg_ = _linearized_bg;
                jacobian.setIdentity();
                covariance.setZero();
                if( VIstaticParm::acc_n == 0 && VIstaticParm::gyr_n == 0 && VIstaticParm::acc_w == 0  && VIstaticParm::gyr_w == 0  )
                    covariance.setIdentity();
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

            Eigen::Matrix<double, 15, 1> evaluate(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi, const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
                                                  const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj, const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj) const
            {
                Eigen::Matrix<double, 15, 1> residuals;

                Eigen::Matrix3d dp_dba = jacobian.block<3, 3>(O_P, O_BA);
                Eigen::Matrix3d dp_dbg = jacobian.block<3, 3>(O_P, O_BG);

                Eigen::Matrix3d dq_dbg = jacobian.block<3, 3>(O_R, O_BG);

                Eigen::Matrix3d dv_dba = jacobian.block<3, 3>(O_V, O_BA);
                Eigen::Matrix3d dv_dbg = jacobian.block<3, 3>(O_V, O_BG);

                Eigen::Vector3d dba = Bai - linearized_ba_;
                Eigen::Vector3d dbg = Bgi - linearized_bg_;

                Eigen::Quaterniond corrected_delta_q = delta_q_ * Utility::deltaQ(dq_dbg * dbg);
                Eigen::Vector3d corrected_delta_v = delta_v_ + dv_dba * dba + dv_dbg * dbg;
                Eigen::Vector3d corrected_delta_p = delta_p_ + dp_dba * dba + dp_dbg * dbg;

                residuals.block<3, 1>(O_P, 0) = Qi.inverse() * (0.5 * VIstaticParm::G_ * sum_dt_ * sum_dt_ + Pj - Pi - Vi * sum_dt_) - corrected_delta_p;
                residuals.block<3, 1>(O_R, 0) = 2 * (corrected_delta_q.inverse() * (Qi.inverse() * Qj)).vec();
                residuals.block<3, 1>(O_V, 0) = Qi.inverse() * (VIstaticParm::G_ * sum_dt_ + Vj - Vi) - corrected_delta_v;
                residuals.block<3, 1>(O_BA, 0) = Baj - Bai;
                residuals.block<3, 1>(O_BG, 0) = Bgj - Bgi;
//                std::cout << "------------IN EVALUATE--------------------" << std::endl;
//                std::cout << "Bai = " << Bai.transpose() << std::endl;
//                std::cout << "Baj = " << Baj.transpose() << std::endl;
//                std::cout << "Bgi = " << Bgi.transpose() << std::endl;
//                std::cout << "Bgj = " << Bgj.transpose() << std::endl;

                return residuals;
            }


            Eigen::Matrix3d GetBg()
            {
                return jacobian.block(3, 12, 3, 3);
            }

            void change_time(const double _t0, const double _t1)
            {
                t0_ = _t0;
                t1_ = _t1;
                sum_dt_ = (static_cast<double>(t1_) - static_cast<double>(t0_)) / 1000.;
                delta_p_.setZero();
                delta_v_.setZero();
                delta_q_.setIdentity();
                jacobian.setIdentity();
                covariance.setZero();
                if( VIstaticParm::acc_n == 0 && VIstaticParm::gyr_n == 0 && VIstaticParm::acc_w == 0  && VIstaticParm::gyr_w == 0  )
                    covariance.setIdentity();
            }

            double sum_dt_; // scond
            double t0_;  // second * 1000
            double t1_;  // second * 1000

            Vec3 delta_p_;
            Vec3 delta_v_;
            Eigen::Quaterniond delta_q_;

            bool good_to_opti_;

            std::vector<double> dt_buf_;
            std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > acc_buf_;
            std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > gyr_buf_;
            Vec3 linearized_ba_, linearized_bg_;
            Vec3 linearized_acc_, linearized_gyr_;

            Vec3 acc_0_, gyr_0_;
            Vec3 acc_1_, gyr_1_;


            Eigen::Matrix<double, 15, 15> jacobian, covariance;
            Eigen::Matrix<double, 18, 18> noise;
        };
    }
}

#endif //OPENMVG_IMU_INTEBASE_HPP
