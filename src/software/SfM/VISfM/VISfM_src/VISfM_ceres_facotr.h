//
// Created by root on 10/16/20.
//

#ifndef OPENMVG_VISFM_CERES_FACOTR_H
#define OPENMVG_VISFM_CERES_FACOTR_H

#include <memory>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"
#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"

//--
//- Define ceres Cost_functor for each OpenMVG camera model
//--

namespace openMVG {
    namespace sfm {

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

    }
}


#endif //OPENMVG_VISFM_CERES_FACOTR_H
