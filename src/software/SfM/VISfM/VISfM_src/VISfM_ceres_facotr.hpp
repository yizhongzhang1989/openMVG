//
// Created by root on 10/16/20.
//

#ifndef OPENMVG_VISFM_CERES_FACOTR_HPP
#define OPENMVG_VISFM_CERES_FACOTR_HPP

#include <memory>
#include <utility>

#include <ceres/ceres.h>
//#include "third_party/ceres-solver/include/ceres/ceres.h"
//#include "third_party/ceres-solver/include/ceres/rotation.h"
#include <ceres/rotation.h>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"
#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"

#include "VISfM_ceres_param.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "Utility.hpp"
#include "VI_static_Parm.hpp"

//--
//- Define ceres Cost_functor for each OpenMVG camera model
//--

namespace openMVG {
    namespace sfm {

        class VISfM_Projection : public ceres::SizedCostFunction<2, 7, 7, 6, 3>
        {
        public:
            VISfM_Projection(Eigen::Vector2d& obs);

            enum : uint8_t {
                OFFSET_FOCAL_LENGTH = 0,
                OFFSET_PRINCIPAL_POINT_X = 1,
                OFFSET_PRINCIPAL_POINT_Y = 2,
                OFFSET_DISTO_K1 = 3,
                OFFSET_DISTO_K2 = 4,
                OFFSET_DISTO_K3 = 5,
            };

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
            void check(double **parameters);
            static Eigen::Vector2d compute_error( double const *parameters_pose, double const *parameters_point, const Eigen::Vector2d& obs, bool& bad_point );
            static Eigen::Matrix2d sqrt_info;
            Eigen::Vector2d point_obs_;
        };

        struct ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3 {
            explicit ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3(const double *const pos_2dpoint)
                    : m_pos_2dpoint(pos_2dpoint) {
            }

            // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
            enum : uint8_t {
                OFFSET_FOCAL_LENGTH = 0,
                OFFSET_PRINCIPAL_POINT_X = 1,
                OFFSET_PRINCIPAL_POINT_Y = 2,
                OFFSET_DISTO_K1 = 3,
                OFFSET_DISTO_K2 = 4,
                OFFSET_DISTO_K3 = 5,
            };

            /**
             * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3 )
             * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
             *   - 3 for rotation(angle axis), 3 for translation
             * @param[in] pos_3dpoint
             * @param[out] out_residuals
             */
            template<typename T>
            bool operator()(
                    const T *const cam_intrinsics,
                    const T *const imu_extrinsics,
                    const T *const cam_imu_extrinsics,
                    const T *const pos_3dpoint,
                    T *out_residuals) const {
                //--
                // Apply external parameters (Pose)
                //--

                const T *Riw = imu_extrinsics;
                Eigen::Map<const Eigen::Matrix<T, 3, 1>> tiw(&imu_extrinsics[3]);
                Eigen::Matrix<T, 3, 1> point_imu;
                ceres::AngleAxisRotatePoint(Riw, pos_3dpoint, point_imu.data());
                point_imu += tiw;

                const T *Rci = cam_imu_extrinsics;
                Eigen::Map<const Eigen::Matrix<T, 3, 1>> tci(&cam_imu_extrinsics[3]);
                Eigen::Matrix<T, 3, 1> point_cam;
                ceres::AngleAxisRotatePoint(Rci, point_imu.data(), point_cam.data());
                point_cam += tci;

                // Transform the point from homogeneous to euclidean (undistorted point)
                const Eigen::Matrix<T, 2, 1> projected_point = point_cam.hnormalized();
                //--
                // Apply intrinsic parameters
                //--

                const T &focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
                const T &principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
                const T &principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
                const T &k1 = cam_intrinsics[OFFSET_DISTO_K1];
                const T &k2 = cam_intrinsics[OFFSET_DISTO_K2];
                const T &k3 = cam_intrinsics[OFFSET_DISTO_K3];

                // Apply distortion (xd,yd) = disto(x_u,y_u)
                const T r2 = projected_point.squaredNorm();
                const T r4 = r2 * r2;
                const T r6 = r4 * r2;
                const T r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);

                Eigen::Map<Eigen::Matrix<T, 2, 1>> residuals(out_residuals);
                residuals << principal_point_x + (projected_point.x() * r_coeff) * focal - m_pos_2dpoint[0],
                        principal_point_y + (projected_point.y() * r_coeff) * focal - m_pos_2dpoint[1];

                return true;
            }

            static int num_residuals() { return 2; }

            // Factory to hide the construction of the CostFunction object from
            // the client code.
            static ceres::CostFunction *Create
                    (
                            const Vec2 &observation,
                            const double weight = 0.0
                    ) {
                if (weight == 0.0) {
                    return
                            (new ceres::AutoDiffCostFunction
                                    <ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 6, 3>(
                                    new ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data())));
                } else {
                    std::runtime_error("not define yet");
                    return
                            (new ceres::AutoDiffCostFunction
                                    <ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 6, 3>(
                                    new ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data())));
//                    return
//                            (new ceres::AutoDiffCostFunction
//                                    <WeightedCostFunction < ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3>, 2, 6, 6,
//                                    3 >
//                                    (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3>
//                                            (new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data()),
//                                             weight)));
                }
            }

            const double *m_pos_2dpoint; // The 2D observation
        };

        class IMUFactor : public ceres::SizedCostFunction<15, 7, 9, 7, 9>
        {
        public:
            IMUFactor() = delete;
//            IMUFactor(std::shared_ptr<IMU_InteBase> _pre_integration)

            IMUFactor(const IMU_InteBase& _pre_integration)
            :pre_integration(_pre_integration)
            {
            }
            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
            {

                Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

                Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
                Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);

                Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

                Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
                Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
                Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);

                Eigen::Map<Eigen::Matrix<double, 15, 1>> residual(residuals);
                residual = pre_integration.evaluate(Pi, Qi, Vi, Bai, Bgi,
                                                     Pj, Qj, Vj, Baj, Bgj);

                Eigen::Matrix<double, 15, 15> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 15, 15>>(pre_integration.covariance.inverse()).matrixL().transpose();
//                Eigen::Matrix<double, 15, 15> sqrt_info;
//                sqrt_info.setIdentity();
                residual = sqrt_info_weight * sqrt_info * residual;
//                std::cout << residual.transpose() << std::endl;

                if (jacobians)
                {
                    double sum_dt = pre_integration.sum_dt_;
                    Eigen::Matrix3d dp_dba = pre_integration.jacobian.template block<3, 3>(O_P, O_BA);
                    Eigen::Matrix3d dp_dbg = pre_integration.jacobian.template block<3, 3>(O_P, O_BG);

                    Eigen::Matrix3d dq_dbg = pre_integration.jacobian.template block<3, 3>(O_R, O_BG);

                    Eigen::Matrix3d dv_dba = pre_integration.jacobian.template block<3, 3>(O_V, O_BA);
                    Eigen::Matrix3d dv_dbg = pre_integration.jacobian.template block<3, 3>(O_V, O_BG);

                    if (jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                        jacobian_pose_i.setZero();

                        jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                        jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * VIstaticParm::G_ * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_ * Utility::deltaQ(dq_dbg * (Bgi - pre_integration.linearized_bg_));
                        jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();

                        jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (VIstaticParm::G_ * sum_dt + Vj - Vi));

                        jacobian_pose_i = sqrt_info * jacobian_pose_i;
                    }
                    if (jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                        jacobian_speedbias_i.setZero();
                        jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;
                        jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
                        jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration.delta_q_).bottomRightCorner<3, 3>() * dq_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();
                        jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
                        jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();

                        jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();

                        jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

                        //ROS_ASSERT(fabs(jacobian_speedbias_i.maxCoeff()) < 1e8);
                        //ROS_ASSERT(fabs(jacobian_speedbias_i.minCoeff()) < 1e8);
                    }
                    if (jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                        jacobian_pose_j.setZero();

                        jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();


                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_ * Utility::deltaQ(dq_dbg * (Bgi - pre_integration.linearized_bg_));
                        jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

                        jacobian_pose_j = sqrt_info * jacobian_pose_j;
                    }
                    if (jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                        jacobian_speedbias_j.setZero();

                        jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

                        jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

                        jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

                        jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;
                    }
                }

                return true;
            }

            //bool Evaluate_Direct(double const *const *parameters, Eigen::Matrix<double, 15, 1> &residuals, Eigen::Matrix<double, 15, 30> &jacobians);

            //void checkCorrection();
            //void checkTransition();
            //void checkJacobian(double **parameters);
            IMU_InteBase pre_integration;
            static Eigen::Matrix<double, 15, 15> sqrt_info_weight;
//            std::shared_ptr<IMU_InteBase> pre_integration;
        };

        class IMUFactorWOBAISE7POSE : public ceres::SizedCostFunction<12, 7, 3, 7, 3>
        {
        public:
            IMUFactorWOBAISE7POSE() = delete;
//            IMUFactor(std::shared_ptr<IMU_InteBase> _pre_integration)

            IMUFactorWOBAISE7POSE(const IMU_InteBase& _pre_integration)
                    :pre_integration(_pre_integration)
            {
            }

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
            {

                Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

                Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);

                Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

                Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);


                Eigen::Vector3d Bai(0.,0.,0.);
                Eigen::Vector3d Bgi(0.,0.,0.);

                Eigen::Vector3d Baj(0.,0.,0.);
                Eigen::Vector3d Bgj(0.,0.,0.);

                Eigen::Map<Eigen::Matrix<double, 12, 1>> residual(residuals);
                residual = pre_integration.evaluate(Pi, Qi, Vi, Bai, Bgi,
                                                            Pj, Qj, Vj, Baj, Bgj).head(12);
//                Eigen::Matrix<double, 12, 12> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 12, 12>>(pre_integration.covariance.inverse()).matrixL().transpose();
                Eigen::Matrix<double, 12, 12> sqrt_info;
                sqrt_info.setIdentity();
                residual = sqrt_info_weight * sqrt_info * residual;
//                std::cout << residual.transpose() << std::endl;

                if (jacobians)
                {
                    double sum_dt = pre_integration.sum_dt_;
                    Eigen::Matrix3d dp_dba = pre_integration.jacobian.template block<3, 3>(O_P, O_BA);
                    Eigen::Matrix3d dp_dbg = pre_integration.jacobian.template block<3, 3>(O_P, O_BG);

                    Eigen::Matrix3d dq_dbg = pre_integration.jacobian.template block<3, 3>(O_R, O_BG);

                    Eigen::Matrix3d dv_dba = pre_integration.jacobian.template block<3, 3>(O_V, O_BA);
                    Eigen::Matrix3d dv_dbg = pre_integration.jacobian.template block<3, 3>(O_V, O_BG);

                    if (jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                        jacobian_pose_i.setZero();

                        jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                        jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * VIstaticParm::G_ * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_ * Utility::deltaQ(dq_dbg * (Bgi - pre_integration.linearized_bg_));
                        jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();

                        jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (VIstaticParm::G_ * sum_dt + Vj - Vi));

                        jacobian_pose_i = sqrt_info * jacobian_pose_i;
                    }
                    if (jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                        jacobian_speedbias_i.setZero();
                        jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;

                        jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();

                        jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;
                    }
                    if (jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                        jacobian_pose_j.setZero();

                        jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();


                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_ * Utility::deltaQ(dq_dbg * (Bgi - pre_integration.linearized_bg_));
                        jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

                        jacobian_pose_j = sqrt_info * jacobian_pose_j;
                    }
                    if (jacobians[5])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                        jacobian_speedbias_j.setZero();

                        jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

                        jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;
                    }
                }

                return true;
            }

            IMU_InteBase pre_integration;
            static Eigen::Matrix<double, 12, 12> sqrt_info_weight;
        };

        class IMUFactorWOBAISE : public ceres::SizedCostFunction<12, 3, 4, 3, 3, 4, 3>
        {
        public:
            IMUFactorWOBAISE() = delete;
//            IMUFactor(std::shared_ptr<IMU_InteBase> _pre_integration)

            IMUFactorWOBAISE(const IMU_InteBase& _pre_integration)
            :pre_integration(_pre_integration)
            {
            }

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
            {

                Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond Qi(parameters[1][3], parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Vector3d Vi(parameters[2][0], parameters[2][1], parameters[2][2]);

                Eigen::Vector3d Pj(parameters[3][0], parameters[3][1], parameters[3][2]);
                Eigen::Quaterniond Qj(parameters[4][3], parameters[4][0], parameters[4][1], parameters[4][2]);
                Eigen::Vector3d Vj(parameters[5][0], parameters[5][1], parameters[5][2]);

                Eigen::Map<Eigen::Matrix<double, 12, 1>> residual(residuals);
                residual = pre_integration.evaluate(Pi, Qi, Vi, Eigen::Vector3d(0.,0.,0.), Eigen::Vector3d(0.,0.,0.),
                                                    Pj, Qj, Vj, Eigen::Vector3d(0.,0.,0.), Eigen::Vector3d(0.,0.,0.)).head(12);

//                Eigen::Matrix<double, 12, 12> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 12, 12>>(pre_integration.covariance.inverse()).matrixL().transpose();
                Eigen::Matrix<double, 12, 12> sqrt_info;
                sqrt_info.setIdentity();
                residual = sqrt_info_weight * sqrt_info * residual;
//                std::cout << residual.transpose() << std::endl;

                if (jacobians)
                {
                    double sum_dt = pre_integration.sum_dt_;
                    Eigen::Matrix3d dp_dba = pre_integration.jacobian.template block<3, 3>(O_P, O_BA);
                    Eigen::Matrix3d dp_dbg = pre_integration.jacobian.template block<3, 3>(O_P, O_BG);

                    Eigen::Matrix3d dq_dbg = pre_integration.jacobian.template block<3, 3>(O_R, O_BG);

                    Eigen::Matrix3d dv_dba = pre_integration.jacobian.template block<3, 3>(O_V, O_BA);
                    Eigen::Matrix3d dv_dbg = pre_integration.jacobian.template block<3, 3>(O_V, O_BG);

                    if (jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_pose_P_i(jacobians[0]);
                        jacobian_pose_P_i.setZero();

                        jacobian_pose_P_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                        jacobian_pose_P_i = sqrt_info * jacobian_pose_P_i;
                    }
                    if (jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 4, Eigen::RowMajor>> jacobian_pose_Q_i(jacobians[1]);
                        jacobian_pose_Q_i.setZero();

                        jacobian_pose_Q_i.block<3, 3>(O_P, O_R - O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * VIstaticParm::G_ * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_;
                        jacobian_pose_Q_i.block<3, 3>(O_R, O_R - O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();

                        jacobian_pose_Q_i.block<3, 3>(O_V, O_R - O_R) = Utility::skewSymmetric(Qi.inverse() * (VIstaticParm::G_ * sum_dt + Vj - Vi));

                        jacobian_pose_Q_i = sqrt_info * jacobian_pose_Q_i;
                    }
                    if (jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[2]);
                        jacobian_speedbias_i.setZero();
                        jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;

                        jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();

                        jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;
                    }
                    if (jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_pose_P_j(jacobians[3]);
                        jacobian_pose_P_j.setZero();

                        jacobian_pose_P_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();

                        jacobian_pose_P_j = sqrt_info * jacobian_pose_P_j;
                    }
                    if (jacobians[4])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 4, Eigen::RowMajor>> jacobian_pose_Q_j(jacobians[4]);
                        jacobian_pose_Q_j.setZero();

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_;
                        jacobian_pose_Q_j.block<3, 3>(O_R, O_R - O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

                        jacobian_pose_Q_j = sqrt_info * jacobian_pose_Q_j;
                    }
                    if (jacobians[5])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[5]);
                        jacobian_speedbias_j.setZero();

                        jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

//                        jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

//                        jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

                        jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;
                    }
                }

                return true;
            }

            /*virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
            {

                Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
                Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);


                Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
                Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);

                Eigen::Map<Eigen::Matrix<double, 12, 1>> residual(residuals);
                residual = pre_integration.evaluate(Pi, Qi, Vi, Eigen::Vector3d(0.,0.,0.), Eigen::Vector3d(0.,0.,0.),
                                                     Pj, Qj, Vj, Eigen::Vector3d(0.,0.,0.), Eigen::Vector3d(0.,0.,0.)).head(12);

//                Eigen::Matrix<double, 12, 12> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 12, 12>>(pre_integration.covariance.inverse()).matrixL().transpose();
                Eigen::Matrix<double, 15, 15> sqrt_info;
                sqrt_info.setIdentity();
                residual = sqrt_info_weight * sqrt_info * residual;
//                std::cout << residual.transpose() << std::endl;

                if (jacobians)
                {
                    double sum_dt = pre_integration.sum_dt_;
                    Eigen::Matrix3d dp_dba = pre_integration.jacobian.template block<3, 3>(O_P, O_BA);
                    Eigen::Matrix3d dp_dbg = pre_integration.jacobian.template block<3, 3>(O_P, O_BG);

                    Eigen::Matrix3d dq_dbg = pre_integration.jacobian.template block<3, 3>(O_R, O_BG);

                    Eigen::Matrix3d dv_dba = pre_integration.jacobian.template block<3, 3>(O_V, O_BA);
                    Eigen::Matrix3d dv_dbg = pre_integration.jacobian.template block<3, 3>(O_V, O_BG);

                    if (jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                        jacobian_pose_i.setZero();

                        jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                        jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * VIstaticParm::G_ * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_;
                        jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();

                        jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (VIstaticParm::G_ * sum_dt + Vj - Vi));

                        jacobian_pose_i = sqrt_info * jacobian_pose_i;
                    }
                    if (jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                        jacobian_speedbias_i.setZero();
                        jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;
//                        jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
//                        jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;

//                        jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration.delta_q_).bottomRightCorner<3, 3>() * dq_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();
//                        jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
//                        jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;

//                        jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();

//                        jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();

                        jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

                        //ROS_ASSERT(fabs(jacobian_speedbias_i.maxCoeff()) < 1e8);
                        //ROS_ASSERT(fabs(jacobian_speedbias_i.minCoeff()) < 1e8);
                    }
                    if (jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                        jacobian_pose_j.setZero();

                        jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();


                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_;
                        jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

                        jacobian_pose_j = sqrt_info * jacobian_pose_j;
                    }
                    if (jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 12, 3, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                        jacobian_speedbias_j.setZero();

                        jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

//                        jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

//                        jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

                        jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;
                    }
                }

                return true;
            }*/

            //bool Evaluate_Direct(double const *const *parameters, Eigen::Matrix<double, 15, 1> &residuals, Eigen::Matrix<double, 15, 30> &jacobians);

            //void checkCorrection();
            //void checkTransition();
            //void checkJacobian(double **parameters);
            IMU_InteBase pre_integration;
            static Eigen::Matrix<double, 12, 12> sqrt_info_weight;
//            std::shared_ptr<IMU_InteBase> pre_integration;
        };


        class IMUFactorBAISE : public ceres::SizedCostFunction<15, 3, 4, 9, 3, 4, 9>
        {
        public:
            IMUFactorBAISE() = delete;
//            IMUFactor(std::shared_ptr<IMU_InteBase> _pre_integration)

            IMUFactorBAISE(const IMU_InteBase& _pre_integration)
                    :pre_integration(_pre_integration)
            {
            }

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
            {

                Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond Qi(parameters[1][3], parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Vector3d Vi(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Vector3d Bai(parameters[2][3], parameters[2][4], parameters[2][5]);
                Eigen::Vector3d Bgi(parameters[2][6], parameters[2][7], parameters[2][8]);


                Eigen::Vector3d Pj(parameters[3][0], parameters[3][1], parameters[3][2]);
                Eigen::Quaterniond Qj(parameters[4][3], parameters[4][0], parameters[4][1], parameters[4][2]);
                Eigen::Vector3d Vj(parameters[5][0], parameters[5][1], parameters[5][2]);
                Eigen::Vector3d Baj(parameters[5][3], parameters[5][4], parameters[5][5]);
                Eigen::Vector3d Bgj(parameters[5][6], parameters[5][7], parameters[5][8]);

                Eigen::Map<Eigen::Matrix<double, 15, 1>> residual(residuals);
                residual = pre_integration.evaluate(Pi, Qi, Vi, Bai, Bgi,
                                                    Pj, Qj, Vj, Baj, Bgj);

//                Eigen::Matrix<double, 12, 12> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 12, 12>>(pre_integration.covariance.inverse()).matrixL().transpose();
                Eigen::Matrix<double, 15, 15> sqrt_info;
                sqrt_info.setIdentity();
                residual = sqrt_info_weight * sqrt_info * residual;
//                std::cout << residual.transpose() << std::endl;

                if (jacobians)
                {
                    double sum_dt = pre_integration.sum_dt_;
                    Eigen::Matrix3d dp_dba = pre_integration.jacobian.template block<3, 3>(O_P, O_BA);
                    Eigen::Matrix3d dp_dbg = pre_integration.jacobian.template block<3, 3>(O_P, O_BG);

                    Eigen::Matrix3d dq_dbg = pre_integration.jacobian.template block<3, 3>(O_R, O_BG);

                    Eigen::Matrix3d dv_dba = pre_integration.jacobian.template block<3, 3>(O_V, O_BA);
                    Eigen::Matrix3d dv_dbg = pre_integration.jacobian.template block<3, 3>(O_V, O_BG);

                    if (jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 3, Eigen::RowMajor>> jacobian_pose_P_i(jacobians[0]);
                        jacobian_pose_P_i.setZero();

                        jacobian_pose_P_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                        jacobian_pose_P_i = sqrt_info * jacobian_pose_P_i;
                    }
                    if (jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 4, Eigen::RowMajor>> jacobian_pose_Q_i(jacobians[1]);
                        jacobian_pose_Q_i.setZero();

                        jacobian_pose_Q_i.block<3, 3>(O_P, O_R - O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * VIstaticParm::G_ * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_;
                        jacobian_pose_Q_i.block<3, 3>(O_R, O_R - O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();

                        jacobian_pose_Q_i.block<3, 3>(O_V, O_R - O_R) = Utility::skewSymmetric(Qi.inverse() * (VIstaticParm::G_ * sum_dt + Vj - Vi));

                        jacobian_pose_Q_i = sqrt_info * jacobian_pose_Q_i;
                    }
                    if (jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[2]);
                        jacobian_speedbias_i.setZero();
                        jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;
                        jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
                        jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration.delta_q_).bottomRightCorner<3, 3>() * dq_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();
                        jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
                        jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;

                        jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();

                        jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();

                        jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;
                    }
                    if (jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 3, Eigen::RowMajor>> jacobian_pose_P_j(jacobians[3]);
                        jacobian_pose_P_j.setZero();

                        jacobian_pose_P_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();

                        jacobian_pose_P_j = sqrt_info * jacobian_pose_P_j;
                    }
                    if (jacobians[4])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 4, Eigen::RowMajor>> jacobian_pose_Q_j(jacobians[4]);
                        jacobian_pose_Q_j.setZero();

                        Eigen::Quaterniond corrected_delta_q = pre_integration.delta_q_;
                        jacobian_pose_Q_j.block<3, 3>(O_R, O_R - O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

                        jacobian_pose_Q_j = sqrt_info * jacobian_pose_Q_j;
                    }
                    if (jacobians[5])
                    {
                        Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[5]);
                        jacobian_speedbias_j.setZero();

                        jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

                        jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

                        jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

                        jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;
                    }
                }

                return true;
            }

            //bool Evaluate_Direct(double const *const *parameters, Eigen::Matrix<double, 15, 1> &residuals, Eigen::Matrix<double, 15, 30> &jacobians);

            //void checkCorrection();
            //void checkTransition();
            //void checkJacobian(double **parameters);
            IMU_InteBase pre_integration;
            static Eigen::Matrix<double, 15, 15> sqrt_info_weight;
//            std::shared_ptr<IMU_InteBase> pre_integration;
        };

        class IMUFactorSoftConstrainsRotation : public ceres::SizedCostFunction<3, 4>
        {
        public:
            IMUFactorSoftConstrainsRotation(const Mat3& _Rwi, double soft );

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
            void check(double **parameters);
            static Eigen::Vector2d compute_error( double const *parameters_pose, double const *parameters_point, const Eigen::Vector2d& obs, bool& bad_point );
//            Eigen::Quaterniond Qwi_;
            std::vector<double> Qwi_double_;
            double soft_;
        };
        class IMUFactorSoftConstrainsTranslation : public ceres::SizedCostFunction<3, 3>
        {

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
            void check(double **parameters);
            static Eigen::Vector2d compute_error( double const *parameters_pose, double const *parameters_point, const Eigen::Vector2d& obs, bool& bad_point );
            Vec3 twi_;
            double soft_;

        public:
            IMUFactorSoftConstrainsTranslation(const Vec3& _twi, double soft );
        };
        class IMUFactorSoftConstrainsVelocity : public ceres::SizedCostFunction<3, 3, 3>
        {
            IMUFactorSoftConstrainsVelocity(const Vec3& _Vi, double soft );

            virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
            void check(double **parameters);
            static Eigen::Vector2d compute_error( double const *parameters_pose, double const *parameters_point, const Eigen::Vector2d& obs, bool& bad_point );
            Vec3 Vi_;
            double soft_;
        };
    }
}


#endif //OPENMVG_VISFM_CERES_FACOTR_HPP
