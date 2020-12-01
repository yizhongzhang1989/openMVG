//
// Created by root on 10/16/20.
//

#include "VISfM_ceres_facotr.hpp"


namespace openMVG
{
    namespace sfm
    {
        Eigen::Matrix<double, 15, 15> IMUFactorBAISE::sqrt_info_weight = Eigen::Matrix<double, 15, 15>::Identity();
        Eigen::Matrix<double, 15, 15> IMUFactor::sqrt_info_weight = Eigen::Matrix<double, 15, 15>::Identity();
        Eigen::Matrix<double, 12, 12> IMUFactorWOBAISE::sqrt_info_weight = Eigen::Matrix<double, 12, 12>::Identity();
        Eigen::Matrix<double, 12, 12> IMUFactorWOBAISE7POSE::sqrt_info_weight = Eigen::Matrix<double, 12, 12>::Identity();
        Eigen::Matrix2d VISfM_Projection::sqrt_info = Eigen::Matrix2d::Identity();
        Eigen::Matrix2d VISfM_ProjectionTd::sqrt_info = Eigen::Matrix2d::Identity();

        VISfM_Projection::VISfM_Projection(Eigen::Vector2d &obs) : point_obs_(obs)
        {

        }

        bool VISfM_Projection::Evaluate(const double *const *parameters, double *residuals, double **jacobians) const
        {
            Eigen::Map<const Eigen::Vector3d> twi(parameters[0]);
            Eigen::Quaterniond Qwi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
            Eigen::Map<const Eigen::Vector3d> tic(parameters[1]);
            Eigen::Quaterniond Qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

            const double &focal = parameters[2][OFFSET_FOCAL_LENGTH];
            const double &principal_point_x = parameters[2][OFFSET_PRINCIPAL_POINT_X];
            const double &principal_point_y = parameters[2][OFFSET_PRINCIPAL_POINT_Y];
            const double &k1 = parameters[2][OFFSET_DISTO_K1];
            const double &k2 = parameters[2][OFFSET_DISTO_K2];
            const double &k3 = parameters[2][OFFSET_DISTO_K3];

            Eigen::Map<const Eigen::Vector3d> pts_w(parameters[3]);

            Eigen::Vector3d pts_imu = Qwi.inverse() * (pts_w - twi);
            Eigen::Vector3d pts_camera = Qic.inverse() * (pts_imu - tic);
            Eigen::Vector2d pts_c_nrom = pts_camera.head(2);
            pts_c_nrom /= pts_camera(2);

            Eigen::Map<Eigen::Vector2d> residual(residuals);
            // Apply distortion (xd,yd) = disto(x_u,y_u)
            const double r2 = pts_c_nrom.squaredNorm();
            const double r4 = r2 * r2;
            const double r6 = r4 * r2;
            const double r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);

            residual << principal_point_x + (pts_c_nrom.x() * r_coeff) * focal - point_obs_(0),
                    principal_point_y + (pts_c_nrom.y() * r_coeff) * focal - point_obs_(1);

            residual = sqrt_info * residual;

//            std::cout << "residuals = " << residual.transpose() << std::endl;

            if( jacobians )
            {
                const double jacobian_rc_ptscnx =
                        2 * k1 * pts_c_nrom.x()
                        + 4 * k2 * r2 * pts_c_nrom.x()
                        + 6 * k3 * r4 * pts_c_nrom.x();
                const double jacobian_rc_ptscny =
                        2 * k1 * pts_c_nrom.y()
                        + 4 * k2 * r2 * pts_c_nrom.y()
                        + 6 * k3 * r4 * pts_c_nrom.y();
                Eigen::Matrix<double, 2, 2> jacobian_e_ptscn;
                jacobian_e_ptscn.setZero();
                jacobian_e_ptscn(0,0) = focal * r_coeff + pts_c_nrom.x() * focal * jacobian_rc_ptscnx;
                jacobian_e_ptscn(0,1) = pts_c_nrom.x() * focal * jacobian_rc_ptscny;
                jacobian_e_ptscn(1,0) = pts_c_nrom.y() * focal * jacobian_rc_ptscnx;
                jacobian_e_ptscn(1,1) = focal * r_coeff + pts_c_nrom.y() * focal * jacobian_rc_ptscny;

                Eigen::Matrix<double, 2, 3> jacobian_ptscn_ptsc;
                jacobian_ptscn_ptsc.setZero();
                jacobian_ptscn_ptsc
                <<
                1. / pts_camera.z(), 0., -1. * pts_camera.x() / ( pts_camera.z() * pts_camera.z() ),
                0., 1. / pts_camera.z(), -1. * pts_camera.y() / ( pts_camera.z() * pts_camera.z() );

                Eigen::Matrix3d Ric = Qic.toRotationMatrix();
                Eigen::Matrix3d Rwi = Qwi.toRotationMatrix();
                Eigen::Matrix3d Identity3d = Eigen::Matrix3d::Identity();

                if(jacobians[0])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 7, Eigen::RowMajor> >jacobian_pose( jacobians[0] );
                    jacobian_pose.setZero();

                    Eigen::Matrix<double, 3, 6> jacobian_ptsc_pose;
                    jacobian_ptsc_pose.block<3,3>(0,0) = - Ric.transpose() * Rwi.transpose();
                    jacobian_ptsc_pose.block<3,3>(0,3) = Ric.transpose() * Utility::skewSymmetric(pts_imu);

                    jacobian_pose.block<2, 6>(0,0) = sqrt_info * jacobian_e_ptscn * jacobian_ptscn_ptsc * jacobian_ptsc_pose;
                }
                if(jacobians[1])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 7, Eigen::RowMajor> >jacobian_ex_pose( jacobians[1] );
                    jacobian_ex_pose.setZero();

                    Eigen::Matrix<double, 3, 6> jacobian_ptsc_ex;
                    jacobian_ptsc_ex.block<3,3>(0,0) = -Ric.transpose();
                    jacobian_ptsc_ex.block<3,3>(0,3) =  Utility::skewSymmetric(pts_camera);

                    jacobian_ex_pose.block<2, 6>(0,0) = sqrt_info * jacobian_e_ptscn * jacobian_ptscn_ptsc * jacobian_ptsc_ex;

//                    jacobian_ex_pose.block<2, 1>(0, 1).setZero();
//                    jacobian_ex_pose.block<2, 1>(0, 2).setZero();

                }
                if(jacobians[2])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 6, Eigen::RowMajor> >jacobian_intrinsics( jacobians[2] );
                    jacobian_intrinsics.setZero();
                    jacobian_intrinsics(0,0) = pts_c_nrom.x() * r_coeff;
                    jacobian_intrinsics(1,0) = pts_c_nrom.y() * r_coeff;
                    jacobian_intrinsics(0,1) = 1;
                    jacobian_intrinsics(1,2) = 1;

                    jacobian_intrinsics(0,3) = pts_c_nrom.x() * focal * r2;
                    jacobian_intrinsics(0,4) = pts_c_nrom.x() * focal * r4;
                    jacobian_intrinsics(0,5) = pts_c_nrom.x() * focal * r6;

                    jacobian_intrinsics(1,3) = pts_c_nrom.y() * focal * r2;
                    jacobian_intrinsics(1,4) = pts_c_nrom.y() * focal * r4;
                    jacobian_intrinsics(1,5) = pts_c_nrom.y() * focal * r6;

                    jacobian_intrinsics = sqrt_info * jacobian_intrinsics;
                }
                if(jacobians[3])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 3, Eigen::RowMajor> >jacobian_point( jacobians[3] );
                    jacobian_point.setZero();

                    Eigen::Matrix3d jacobian_ptsc_ptsw = Ric.transpose() * Rwi.transpose();

                    jacobian_point = sqrt_info * jacobian_e_ptscn * jacobian_ptscn_ptsc * jacobian_ptsc_ptsw;
                }
                bool numb_jaco = false;
                if( numb_jaco )
                {
                    const double eps = 1e-6;
                    Eigen::Matrix<double, 2, 21> num_jacobian;
                    for (int k = 0; k < 21; k++)
                    {
                        Eigen::Vector3d twi_numb( parameters[0][0], parameters[0][1], parameters[0][2] );
                        Eigen::Quaterniond Qwi_numb(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
                        Eigen::Vector3d tic_numb( parameters[1][0], parameters[1][1], parameters[1][2] );
                        Eigen::Quaterniond Qic_numb(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                        double focal_numb = parameters[2][OFFSET_FOCAL_LENGTH];
                        double principal_point_x_numb = parameters[2][OFFSET_PRINCIPAL_POINT_X];
                        double principal_point_y_numb = parameters[2][OFFSET_PRINCIPAL_POINT_Y];
                        double k1_numb = parameters[2][OFFSET_DISTO_K1];
                        double k2_numb = parameters[2][OFFSET_DISTO_K2];
                        double k3_numb = parameters[2][OFFSET_DISTO_K3];

                        Eigen::Vector3d pts_w_numb( parameters[3][0], parameters[3][1], parameters[3][2] );


                        int a = k / 3, b = k % 3;
                        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

                        if (a == 0)
                            twi_numb += delta;
                        else if (a == 1)
                            Qwi_numb = Qwi_numb * Utility::deltaQ(delta);
                        else if (a == 2)
                            tic_numb += delta;
                        else if (a == 3)
                            Qic_numb = Qic_numb * Utility::deltaQ(delta);
                        else if (a == 4)
                        {
                            switch (b) {
                                case 0:
                                    focal_numb += eps;
                                    break;
                                case 1:
                                    principal_point_x_numb += eps;
                                    break;
                                case 2:
                                    principal_point_y_numb += eps;
                                    break;
                                default:
                                    break;
                            }
                        }
                        else if(a==5)
                        {
                            switch (b) {
                                case 0:
                                    k1_numb += eps;
                                    break;
                                case 1:
                                    k2_numb += eps;
                                    break;
                                case 2:
                                    k3_numb += eps;
                                    break;
                                default:
                                    break;
                            }
                        }
                        else if (a == 6)
                            pts_w_numb += delta;


                        Eigen::Vector3d pts_imu_numb = Qwi_numb.inverse() * (pts_w_numb - twi_numb);
                        Eigen::Vector3d pts_camera_numb = Qic_numb.inverse() * (pts_imu_numb - tic_numb);
                        Eigen::Vector2d pts_c_nrom_numb = pts_camera_numb.head(2);
                        pts_c_nrom_numb /= pts_camera_numb(2);

                        // Apply distortion (xd,yd) = disto(x_u,y_u)
                        const double r2_numb = pts_c_nrom_numb.squaredNorm();
                        const double r4_numb = r2_numb * r2_numb;
                        const double r6_numb = r4_numb * r2_numb;
                        const double r_coeff_numb = (1.0 + k1_numb * r2_numb + k2_numb * r4_numb + k3_numb * r6_numb);

                        Eigen::Vector2d residual_numb;
                        residual_numb << principal_point_x_numb + (pts_c_nrom_numb.x() * r_coeff_numb) * focal_numb - point_obs_(0),
                                principal_point_y_numb + (pts_c_nrom_numb.y() * r_coeff_numb) * focal_numb - point_obs_(1);

                        residual_numb = sqrt_info * residual_numb;
                        num_jacobian.col(k) = (residual_numb - residual) / eps;
                    }

//                    std::cout << "========================" << std::endl;
                    if(jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
                        jacobian_pose.setZero();
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_pose " << std::endl << jacobian_pose.block<2, 6>(0,0) << std::endl;
//                        std::cout << point_obs_.transpose() << " jacobian_pose num " << num_jacobian.block<2, 6>( 0, 0 ) << std::endl;
                        jacobian_pose.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 0 );
                    }
                    if(jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex(jacobians[1]);
                        jacobian_ex.setZero();
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_ex " << std::endl <<   jacobian_ex.block<2, 6>(0, 0) << std::endl;
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_ex num " << std::endl <<  num_jacobian.block<2, 6>( 0, 6 ) << std::endl;
                        jacobian_ex.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 6 );
                    }
                    if(jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> jacobian_intrinsics(jacobians[2]);
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_intrinsics " << std::endl <<  jacobian_intrinsics << std::endl;
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_intrinsics num " << std::endl <<  num_jacobian.block<2, 6>( 0, 12 ) << std::endl;
                        jacobian_intrinsics.setZero();
                        jacobian_intrinsics.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 12 );
                    }
                    if(jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_point(jacobians[3]);
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_point " << std::endl << jacobian_point << std::endl;
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_point num " << std::endl << num_jacobian.block<2, 3>( 0, 18 ) << std::endl;
                        jacobian_point.setZero();
                        jacobian_point.block<2, 3>(0, 0) = num_jacobian.block<2, 3>( 0, 18 );
                    }
//                    std::cout << "========================" << std::endl;
                }
            }

            return true;
        }

        VISfM_ProjectionTd::VISfM_ProjectionTd(const Eigen::Vector2d& obs,
                                               const Eigen::Vector2d& point_velocity)
                                               : point_obs_(obs), point_velocity_(point_velocity)
        {

        }

        bool VISfM_ProjectionTd::Evaluate(const double *const *parameters, double *residuals, double **jacobians) const
        {
            Eigen::Map<const Eigen::Vector3d> twi(parameters[0]);
            Eigen::Quaterniond Qwi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
            Eigen::Map<const Eigen::Vector3d> tic(parameters[1]);
            Eigen::Quaterniond Qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

            const double &focal = parameters[2][OFFSET_FOCAL_LENGTH];
            const double &principal_point_x = parameters[2][OFFSET_PRINCIPAL_POINT_X];
            const double &principal_point_y = parameters[2][OFFSET_PRINCIPAL_POINT_Y];
            const double &k1 = parameters[2][OFFSET_DISTO_K1];
            const double &k2 = parameters[2][OFFSET_DISTO_K2];
            const double &k3 = parameters[2][OFFSET_DISTO_K3];

            Eigen::Map<const Eigen::Vector3d> pts_w(parameters[3]);

            double td = parameters[4][0];

            Eigen::Vector2d point_obs_td = point_obs_ - td * point_velocity_;

            Eigen::Vector3d pts_imu = Qwi.inverse() * (pts_w - twi);
            Eigen::Vector3d pts_camera = Qic.inverse() * (pts_imu - tic);
            Eigen::Vector2d pts_c_nrom = pts_camera.head(2);
            pts_c_nrom /= pts_camera(2);

            Eigen::Map<Eigen::Vector2d> residual(residuals);
            // Apply distortion (xd,yd) = disto(x_u,y_u)
            const double r2 = pts_c_nrom.squaredNorm();
            const double r4 = r2 * r2;
            const double r6 = r4 * r2;
            const double r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);

            residual << principal_point_x + (pts_c_nrom.x() * r_coeff) * focal - point_obs_td(0),
                    principal_point_y + (pts_c_nrom.y() * r_coeff) * focal - point_obs_td(1);

            residual = sqrt_info * residual;

//            std::cout << "residuals = " << residual.transpose() << std::endl;

            if( jacobians )
            {
                const double jacobian_rc_ptscnx =
                        2 * k1 * pts_c_nrom.x()
                        + 4 * k2 * r2 * pts_c_nrom.x()
                        + 6 * k3 * r4 * pts_c_nrom.x();
                const double jacobian_rc_ptscny =
                        2 * k1 * pts_c_nrom.y()
                        + 4 * k2 * r2 * pts_c_nrom.y()
                        + 6 * k3 * r4 * pts_c_nrom.y();
                Eigen::Matrix<double, 2, 2> jacobian_e_ptscn;
                jacobian_e_ptscn.setZero();
                jacobian_e_ptscn(0,0) = focal * r_coeff + pts_c_nrom.x() * focal * jacobian_rc_ptscnx;
                jacobian_e_ptscn(0,1) = pts_c_nrom.x() * focal * jacobian_rc_ptscny;
                jacobian_e_ptscn(1,0) = pts_c_nrom.y() * focal * jacobian_rc_ptscnx;
                jacobian_e_ptscn(1,1) = focal * r_coeff + pts_c_nrom.y() * focal * jacobian_rc_ptscny;

                Eigen::Matrix<double, 2, 3> jacobian_ptscn_ptsc;
                jacobian_ptscn_ptsc.setZero();
                jacobian_ptscn_ptsc
                        <<
                        1. / pts_camera.z(), 0., -1. * pts_camera.x() / ( pts_camera.z() * pts_camera.z() ),
                        0., 1. / pts_camera.z(), -1. * pts_camera.y() / ( pts_camera.z() * pts_camera.z() );

                Eigen::Matrix3d Ric = Qic.toRotationMatrix();
                Eigen::Matrix3d Rwi = Qwi.toRotationMatrix();
                Eigen::Matrix3d Identity3d = Eigen::Matrix3d::Identity();

                if(jacobians[0])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 7, Eigen::RowMajor> >jacobian_pose( jacobians[0] );
                    jacobian_pose.setZero();

                    Eigen::Matrix<double, 3, 6> jacobian_ptsc_pose;
                    jacobian_ptsc_pose.block<3,3>(0,0) = - Ric.transpose() * Rwi.transpose();
                    jacobian_ptsc_pose.block<3,3>(0,3) = Ric.transpose() * Utility::skewSymmetric(pts_imu);

                    jacobian_pose.block<2, 6>(0,0) = sqrt_info * jacobian_e_ptscn * jacobian_ptscn_ptsc * jacobian_ptsc_pose;
                }
                if(jacobians[1])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 7, Eigen::RowMajor> >jacobian_ex_pose( jacobians[1] );
                    jacobian_ex_pose.setZero();

                    Eigen::Matrix<double, 3, 6> jacobian_ptsc_ex;
                    jacobian_ptsc_ex.block<3,3>(0,0) = -Ric.transpose();
                    jacobian_ptsc_ex.block<3,3>(0,3) =  Utility::skewSymmetric(pts_camera);

                    jacobian_ex_pose.block<2, 6>(0,0) = sqrt_info * jacobian_e_ptscn * jacobian_ptscn_ptsc * jacobian_ptsc_ex;

//                    jacobian_ex_pose.block<2, 1>(0, 1).setZero();
//                    jacobian_ex_pose.block<2, 1>(0, 2).setZero();

                }
                if(jacobians[2])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 6, Eigen::RowMajor> >jacobian_intrinsics( jacobians[2] );
                    jacobian_intrinsics.setZero();
                    jacobian_intrinsics(0,0) = pts_c_nrom.x() * r_coeff;
                    jacobian_intrinsics(1,0) = pts_c_nrom.y() * r_coeff;
                    jacobian_intrinsics(0,1) = 1;
                    jacobian_intrinsics(1,2) = 1;

                    jacobian_intrinsics(0,3) = pts_c_nrom.x() * focal * r2;
                    jacobian_intrinsics(0,4) = pts_c_nrom.x() * focal * r4;
                    jacobian_intrinsics(0,5) = pts_c_nrom.x() * focal * r6;

                    jacobian_intrinsics(1,3) = pts_c_nrom.y() * focal * r2;
                    jacobian_intrinsics(1,4) = pts_c_nrom.y() * focal * r4;
                    jacobian_intrinsics(1,5) = pts_c_nrom.y() * focal * r6;

                    jacobian_intrinsics = sqrt_info * jacobian_intrinsics;
                }
                if(jacobians[3])
                {
                    Eigen::Map< Eigen::Matrix<double, 2, 3, Eigen::RowMajor> >jacobian_point( jacobians[3] );
                    jacobian_point.setZero();

                    Eigen::Matrix3d jacobian_ptsc_ptsw = Ric.transpose() * Rwi.transpose();

                    jacobian_point = sqrt_info * jacobian_e_ptscn * jacobian_ptscn_ptsc * jacobian_ptsc_ptsw;
                }
                if(jacobians[4])
                {
                    Eigen::Map<Eigen::Vector2d> jacobian_td(jacobians[4]);
                    jacobian_td.setZero();

                    jacobian_td  = sqrt_info * point_velocity_;
                }
                bool numb_jaco = false;
                if( numb_jaco )
                {
                    const double eps = 1e-6;
                    Eigen::Matrix<double, 2, 21> num_jacobian;
                    for (int k = 0; k < 21; k++)
                    {
                        Eigen::Vector3d twi_numb( parameters[0][0], parameters[0][1], parameters[0][2] );
                        Eigen::Quaterniond Qwi_numb(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
                        Eigen::Vector3d tic_numb( parameters[1][0], parameters[1][1], parameters[1][2] );
                        Eigen::Quaterniond Qic_numb(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                        double focal_numb = parameters[2][OFFSET_FOCAL_LENGTH];
                        double principal_point_x_numb = parameters[2][OFFSET_PRINCIPAL_POINT_X];
                        double principal_point_y_numb = parameters[2][OFFSET_PRINCIPAL_POINT_Y];
                        double k1_numb = parameters[2][OFFSET_DISTO_K1];
                        double k2_numb = parameters[2][OFFSET_DISTO_K2];
                        double k3_numb = parameters[2][OFFSET_DISTO_K3];

                        Eigen::Vector3d pts_w_numb( parameters[3][0], parameters[3][1], parameters[3][2] );


                        int a = k / 3, b = k % 3;
                        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

                        if (a == 0)
                            twi_numb += delta;
                        else if (a == 1)
                            Qwi_numb = Qwi_numb * Utility::deltaQ(delta);
                        else if (a == 2)
                            tic_numb += delta;
                        else if (a == 3)
                            Qic_numb = Qic_numb * Utility::deltaQ(delta);
                        else if (a == 4)
                        {
                            switch (b) {
                                case 0:
                                    focal_numb += eps;
                                    break;
                                case 1:
                                    principal_point_x_numb += eps;
                                    break;
                                case 2:
                                    principal_point_y_numb += eps;
                                    break;
                                default:
                                    break;
                            }
                        }
                        else if(a==5)
                        {
                            switch (b) {
                                case 0:
                                    k1_numb += eps;
                                    break;
                                case 1:
                                    k2_numb += eps;
                                    break;
                                case 2:
                                    k3_numb += eps;
                                    break;
                                default:
                                    break;
                            }
                        }
                        else if (a == 6)
                            pts_w_numb += delta;


                        Eigen::Vector3d pts_imu_numb = Qwi_numb.inverse() * (pts_w_numb - twi_numb);
                        Eigen::Vector3d pts_camera_numb = Qic_numb.inverse() * (pts_imu_numb - tic_numb);
                        Eigen::Vector2d pts_c_nrom_numb = pts_camera_numb.head(2);
                        pts_c_nrom_numb /= pts_camera_numb(2);

                        // Apply distortion (xd,yd) = disto(x_u,y_u)
                        const double r2_numb = pts_c_nrom_numb.squaredNorm();
                        const double r4_numb = r2_numb * r2_numb;
                        const double r6_numb = r4_numb * r2_numb;
                        const double r_coeff_numb = (1.0 + k1_numb * r2_numb + k2_numb * r4_numb + k3_numb * r6_numb);

                        Eigen::Vector2d residual_numb;
                        residual_numb << principal_point_x_numb + (pts_c_nrom_numb.x() * r_coeff_numb) * focal_numb - point_obs_(0),
                                principal_point_y_numb + (pts_c_nrom_numb.y() * r_coeff_numb) * focal_numb - point_obs_(1);

                        residual_numb = sqrt_info * residual_numb;
                        num_jacobian.col(k) = (residual_numb - residual) / eps;
                    }

//                    std::cout << "========================" << std::endl;
                    if(jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
                        jacobian_pose.setZero();
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_pose " << std::endl << jacobian_pose.block<2, 6>(0,0) << std::endl;
//                        std::cout << point_obs_.transpose() << " jacobian_pose num " << num_jacobian.block<2, 6>( 0, 0 ) << std::endl;
                        jacobian_pose.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 0 );
                    }
                    if(jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex(jacobians[1]);
                        jacobian_ex.setZero();
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_ex " << std::endl <<   jacobian_ex.block<2, 6>(0, 0) << std::endl;
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_ex num " << std::endl <<  num_jacobian.block<2, 6>( 0, 6 ) << std::endl;
                        jacobian_ex.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 6 );
                    }
                    if(jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> jacobian_intrinsics(jacobians[2]);
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_intrinsics " << std::endl <<  jacobian_intrinsics << std::endl;
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_intrinsics num " << std::endl <<  num_jacobian.block<2, 6>( 0, 12 ) << std::endl;
                        jacobian_intrinsics.setZero();
                        jacobian_intrinsics.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 12 );
                    }
                    if(jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_point(jacobians[3]);
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_point " << std::endl << jacobian_point << std::endl;
//                        std::cout << point_obs_.transpose() << std::endl << " jacobian_point num " << std::endl << num_jacobian.block<2, 3>( 0, 18 ) << std::endl;
                        jacobian_point.setZero();
                        jacobian_point.block<2, 3>(0, 0) = num_jacobian.block<2, 3>( 0, 18 );
                    }
//                    std::cout << "========================" << std::endl;
                }
            }

            return true;
        }

        TdRegularizationTerm::TdRegularizationTerm()
        {

        }

        bool TdRegularizationTerm::Evaluate(const double *const *parameters, double *residuals, double **jacobians) const
        {
            double td_i = parameters[0][0];
            double td_j = parameters[1][0];

            *residuals = td_i - td_j;

            if( jacobians )
            {

                if(jacobians[0])
                {
                    jacobians[0][0] = 1;
                }
                if(jacobians[1])
                {

                    jacobians[1][0] = -1;

                }
            }

            return true;
        }


        IMUFactorSoftConstrainsTranslation::IMUFactorSoftConstrainsTranslation(const Vec3 &_twi, double soft)
        {
            twi_ = _twi;
            soft_ = soft;
        }

        bool IMUFactorSoftConstrainsTranslation::Evaluate(const double *const *parameters, double *residuals, double **jacobians) const
        {
            Eigen::Vector3d twi(parameters[0][0], parameters[0][1], parameters[0][2]);

            Eigen::Map<Eigen::Vector3d> residual(residuals);
            residual = twi - twi_;
            residual = soft_ * residual;

            if( jacobians )
            {
                if(jacobians[0])
                {
                    Eigen::Map< Eigen::Matrix<double, 3, 3, Eigen::RowMajor> >jacobian_translation( jacobians[0] );
                    jacobian_translation.setIdentity();
                    jacobian_translation *= soft_;
                }

                bool numb_jaco = false;
                if( numb_jaco )
                {
                    const double eps = 1e-6;
                    Eigen::Matrix<double, 3, 3> num_jacobian;
                    for (int k = 0; k < 3; k++)
                    {
                        Eigen::Vector3d twi_numb( parameters[0][0], parameters[0][1], parameters[0][2] );

                        int a = k / 3, b = k % 3;
                        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

                        if (a == 0)
                            twi_numb += delta;

                        Eigen::Vector3d residual_numb;
                        residual_numb = twi_numb - twi_;
                        residual_numb = soft_ * residual_numb;
                        num_jacobian.col(k) = (residual_numb - residual) / eps;
                    }
                    std::cout << "========================" << std::endl;
                    if(jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> jacobian_translation(jacobians[0]);

                        std::cout << jacobian_translation << std::endl;
                        std::cout << "------------------------" << std::endl;
                        std::cout << num_jacobian << std::endl;
                    }
                    std::cout << "========================" << std::endl;
                }
            }

            return true;
        }

        IMUFactorSoftConstrainsRotation::IMUFactorSoftConstrainsRotation(const Mat3 &_Rwi, double soft)
        {
            Eigen::Quaterniond Qwi(_Rwi);
            Qwi_double_ = { Qwi.w(), Qwi.x(), Qwi.y(), Qwi.z() };
            soft_ = soft;
        }

        bool IMUFactorSoftConstrainsRotation::Evaluate(const double *const *parameters, double *residuals, double **jacobians) const
        {
            Eigen::Quaterniond Qwi_(Qwi_double_[0], Qwi_double_[1], Qwi_double_[2], Qwi_double_[3]);
            Eigen::Quaterniond Qwi(parameters[0][3], parameters[0][0], parameters[0][1], parameters[0][2]);

            Eigen::Map<Eigen::Vector3d> residual(residuals);

            residual = 2 * (Qwi_.inverse() * Qwi).vec();
            residual = soft_ * residual;

            if( jacobians )
            {
                if(jacobians[0])
                {
                    Eigen::Map< Eigen::Matrix<double, 3, 4, Eigen::RowMajor> >jacobian_rotation( jacobians[0] );
                    jacobian_rotation.setZero();
                    jacobian_rotation.block<3,3>(0,0) = Utility::Qleft(Qwi_.inverse() * Qwi).bottomRightCorner<3, 3>();
                    jacobian_rotation *= soft_;
                }

                bool numb_jaco = false;
                if( numb_jaco )
                {
                    const double eps = 1e-6;
                    Eigen::Matrix<double, 3, 3> num_jacobian;
                    for (int k = 0; k < 3; k++)
                    {
                        Eigen::Quaterniond Qwi_numb(parameters[0][3], parameters[0][0], parameters[0][1], parameters[0][2]);
                        int a = k / 3, b = k % 3;
                        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                        if (a == 0)
                            Qwi_numb = Qwi_numb * Utility::deltaQ(delta);

                        Eigen::Vector3d residual_numb;
                        residual_numb = 2 * (Qwi_.inverse() * Qwi_numb).vec();
                        residual_numb = soft_ * residual_numb;
                        num_jacobian.col(k) = (residual_numb - residual) / eps;
                    }

                    std::cout << "========================" << std::endl;
                    if(jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> jacobian_rotation(jacobians[0]);

                        std::cout << jacobian_rotation << std::endl;
                        std::cout << "------------------------" << std::endl;
                        std::cout << num_jacobian << std::endl;
                    }
                    std::cout << "========================" << std::endl;
                }
            }

            return true;


        }

    }
}