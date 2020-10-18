//
// Created by root on 10/16/20.
//

#include "VISfM_ceres_facotr.hpp"


namespace openMVG
{
    namespace sfm
    {
        Eigen::Matrix2d VISfM_Projection::sqrt_info = Eigen::Matrix2d::Identity();

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

            if( jacobians )
            {
                if(jacobians[0])
                {

                }
                if(jacobians[1])
                {

                }
                if(jacobians[2])
                {

                }
                if(jacobians[3])
                {

                }
                bool numb_jaco = true;
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
                    if(jacobians[0])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
                        jacobian_pose.setZero();
                        jacobian_pose.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 0 );
                    }
                    if(jacobians[1])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex(jacobians[1]);
                        jacobian_ex.setZero();
                        jacobian_ex.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 6 );
                    }
                    if(jacobians[2])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> jacobian_intrinsics(jacobians[2]);
                        jacobian_intrinsics.setZero();
                        jacobian_intrinsics.block<2, 6>(0, 0) = num_jacobian.block<2, 6>( 0, 12 );
                    }
                    if(jacobians[3])
                    {
                        Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_point(jacobians[3]);
                        jacobian_point.setZero();
                        jacobian_point.block<2, 3>(0, 0) = num_jacobian.block<2, 3>( 0, 18 );
                    }
                }
            }

            return true;
        }



    }
}