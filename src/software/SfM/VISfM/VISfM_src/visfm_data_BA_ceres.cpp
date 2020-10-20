//
// Created by root on 10/15/20.
//

#include "visfm_data_BA_ceres.hpp"


#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include "ceres/problem.h"
#include "ceres/solver.h"
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
//- Robust estimation - LMeds (since no threshold can be defined)
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include <ceres/rotation.h>
#include <ceres/types.h>

#include <iostream>
#include <limits>


namespace openMVG
{
    namespace sfm
    {
        using namespace openMVG::cameras;
        using namespace openMVG::geometry;

        // Ceres CostFunctor used for SfM pose center to GPS pose center minimization
        struct PoseCenterConstraintCostFunction
        {
            Vec3 weight_;
            Vec3 pose_center_constraint_;

            PoseCenterConstraintCostFunction
                    (
                            const Vec3 & center,
                            const Vec3 & weight
                    ): weight_(weight), pose_center_constraint_(center)
            {
            }

            template <typename T> bool
            operator()
                    (
                            const T* const cam_extrinsics, // R_t
                            T* residuals
                    )
            const
            {
                using Vec3T = Eigen::Matrix<T,3,1>;
                Eigen::Map<const Vec3T> cam_R(&cam_extrinsics[0]);
                Eigen::Map<const Vec3T> cam_t(&cam_extrinsics[3]);
                const Vec3T cam_R_transpose(-cam_R);

                Vec3T pose_center;
                // Rotate the point according the camera rotation
                ceres::AngleAxisRotatePoint(cam_R_transpose.data(), cam_t.data(), pose_center.data());
                pose_center = pose_center * T(-1);

                Eigen::Map<Vec3T> residuals_eigen(residuals);
                residuals_eigen = weight_.cast<T>().cwiseProduct(pose_center - pose_center_constraint_.cast<T>());

                return true;
            }
        };


        Bundle_Adjustment_IMU_Ceres::BA_Ceres_options::BA_Ceres_options
                (
                        const bool bVerbose,
                        bool bmultithreaded
                )
                : bVerbose_(bVerbose),
                  nb_threads_(1),
                  parameter_tolerance_(1e-8), //~= numeric_limits<float>::epsilon()
                  bUse_loss_function_(true)
        {
#ifdef OPENMVG_USE_OPENMP
            nb_threads_ = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP
            if (!bmultithreaded)
                nb_threads_ = 1;

            bCeres_summary_ = false;

            // Default configuration use a DENSE representation
            linear_solver_type_ = ceres::DENSE_SCHUR;
            preconditioner_type_ = ceres::JACOBI;
            // If Sparse linear solver are available
            // Descending priority order by efficiency (SUITE_SPARSE > CX_SPARSE > EIGEN_SPARSE)
            if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
            {
                sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
                linear_solver_type_ = ceres::SPARSE_SCHUR;
            }
            else
            {
                if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
                {
                    sparse_linear_algebra_library_type_ = ceres::CX_SPARSE;
                    linear_solver_type_ = ceres::SPARSE_SCHUR;
                }
                else
                if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
                {
                    sparse_linear_algebra_library_type_ = ceres::EIGEN_SPARSE;
                    linear_solver_type_ = ceres::SPARSE_SCHUR;
                }
            }
        }

        Bundle_Adjustment_IMU_Ceres::Bundle_Adjustment_IMU_Ceres
                (
                        const Bundle_Adjustment_IMU_Ceres::BA_Ceres_options & options
                )
                : ceres_options_(options)
        {}

        Bundle_Adjustment_IMU_Ceres::BA_Ceres_options &
        Bundle_Adjustment_IMU_Ceres::ceres_options()
        {
            return ceres_options_;
        }

        bool Bundle_Adjustment_IMU_Ceres::Adjust_onlyvisual(sfm::SfM_Data &sfm_data, const Optimize_Options &options)
        {
            if( options.use_motion_priors_opt )
            {
                std::cerr << "not define yet" << std::endl;
                assert(0);
            }

            ceres::Problem problem;

            double ex_paparm[7];
            {
                Mat3 Ric = sfm_data.IG_Ric;
                Eigen::Quaterniond Qic(Ric);
                Vec3 tic = sfm_data.IG_tic;
                ex_paparm[0] = tic(0);
                ex_paparm[1] = tic(1);
                ex_paparm[2] = tic(2);
                ex_paparm[3] = Qic.x();
                ex_paparm[4] = Qic.y();
                ex_paparm[5] = Qic.z();
                ex_paparm[6] = Qic.w();

                ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
                problem.AddParameterBlock(ex_paparm, 7, local_parameterization);  // p,q
                problem.SetParameterBlockConstant(ex_paparm);
            }

            // Data wrapper for refinement:
            Hash_Map<IndexT, std::vector<double>> map_intrinsics;
            Hash_Map<IndexT, std::vector<double>> map_poses;

            // Setup Poses data & subparametrization
            for (const auto & pose_it : sfm_data.poses)
            {
                const IndexT indexPose = pose_it.first;

                const Pose3 & pose = pose_it.second;
                const Mat3 Rwc = pose.rotation().transpose();
                const Vec3 twc = pose.center();

                Vec3 tci = - sfm_data.IG_Ric.transpose() * sfm_data.IG_tic;
                Mat3 Rwi = Rwc * sfm_data.IG_Ric.transpose();
                Eigen::Quaterniond Qwi(Rwi);
                Vec3 twi = twc + Rwc * tci;

                // angleAxis + translation
                map_poses[indexPose] = {
                        twi(0),
                        twi(1),
                        twi(2),
                        Qwi.x(),
                        Qwi.y(),
                        Qwi.z(),
                        Qwi.w()
                };

                double * parameter_block = &map_poses.at(indexPose)[0];
                ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
                problem.AddParameterBlock(parameter_block, 7, local_parameterization);  // p,q
//                problem.AddParameterBlock(parameter_block, 6);
                if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
                {
                    // set the whole parameter block as constant for best performance
                    problem.SetParameterBlockConstant(parameter_block);
                }
                else  // Subset parametrization
                {
                    std::vector<int> vec_constant_extrinsic;
                    // If we adjust only the translation, we must set ROTATION as constant
                    if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
                    {
                        // Subset rotation parametrization
                        vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), {0,1,2});
                    }
                    // If we adjust only the rotation, we must set TRANSLATION as constant
                    if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
                    {
                        // Subset translation parametrization
                        vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), {3,4,5,6});
                    }
                    if (!vec_constant_extrinsic.empty())
                    {
                        assert(0);
                        ceres::SubsetParameterization *subset_parameterization =
                                new ceres::SubsetParameterization(7, vec_constant_extrinsic);
//                        auto sssssssson =
//                                new ceres::SubsetParameterization(7, vec_constant_extrinsic, local_parameterization);
                        problem.SetParameterization(parameter_block, subset_parameterization);
                    }
                }
            }

            // Setup Intrinsics data & subparametrization
            for (const auto & intrinsic_it : sfm_data.intrinsics)
            {
                const IndexT indexCam = intrinsic_it.first;

                if (isValid(intrinsic_it.second->getType()))
                {
                    map_intrinsics[indexCam] = intrinsic_it.second->getParams();
                    if (!map_intrinsics.at(indexCam).empty())
                    {
                        double * parameter_block = &map_intrinsics.at(indexCam)[0];
                        problem.AddParameterBlock(parameter_block, map_intrinsics.at(indexCam).size());
                        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
                        {
                            // set the whole parameter block as constant for best performance
                            problem.SetParameterBlockConstant(parameter_block);
                        }
                        else
                        {
                            const std::vector<int> vec_constant_intrinsic =
                                    intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
                            if (!vec_constant_intrinsic.empty())
                            {
                                ceres::SubsetParameterization *subset_parameterization =
                                        new ceres::SubsetParameterization(
                                                map_intrinsics.at(indexCam).size(), vec_constant_intrinsic);
                                problem.SetParameterization(parameter_block, subset_parameterization);
                            }
                        }
                    }
                }
                else
                {
                    std::cerr << "Unsupported camera type." << std::endl;
                }
            }

            // Set a LossFunction to be less penalized by false measurements
            //  - set it to nullptr if you don't want use a lossFunction.
            ceres::LossFunction * p_LossFunction =
                    ceres_options_.bUse_loss_function_ ?
                    new ceres::HuberLoss(Square(4.0))
                                                       : nullptr;

            // For all visibility add reprojections errors:
            for (auto & structure_landmark_it : sfm_data.structure)
            {
                const Observations & obs = structure_landmark_it.second.obs;

                for (const auto & obs_it : obs)
                {
                    // Build the residual block corresponding to the track observation:
                    const View * view = sfm_data.views.at(obs_it.first).get();

                    // Each Residual block takes a point and a camera as input and outputs a 2
                    // dimensional residual. Internally, the cost function stores the observed
                    // image location and compares the reprojection against the observation.
//                    ceres::CostFunction* cost_function =
//                            (new ceres::AutoDiffCostFunction
//                                    <ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 6, 3>(
//                                    new ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3(obs_it.second.x.data())));
//                            assert( sfm_data.intrinsics.at(view->id_intrinsic).get()->getType() == PINHOLE_CAMERA_RADIAL3 );
//                            IntrinsicsToCostFunction(sfm_data.intrinsics.at(view->id_intrinsic).get(),
//                                                     obs_it.second.x);
                    Eigen::Vector2d ob_i = obs_it.second.x;
                    VISfM_Projection* cost_function = new VISfM_Projection( ob_i );

                    if (cost_function)
                    {
                        if (!map_intrinsics.at(view->id_intrinsic).empty())
                        {
                            problem.AddResidualBlock(cost_function,
                                                     p_LossFunction,
                                                     &map_poses.at(view->id_pose)[0],
                                                     ex_paparm,
                                                     &map_intrinsics.at(view->id_intrinsic)[0],
                                                     structure_landmark_it.second.X.data());

//                            problem.AddResidualBlock(cost_function,
//                                                     p_LossFunction,
//                                                     &map_intrinsics.at(view->id_intrinsic)[0],
//                                                     &map_poses.at(view->id_pose)[0],
//                                                     ex_paparm,
//                                                     structure_landmark_it.second.X.data());
                        }
                        else
                        {
                            assert(0);
//                            problem.AddResidualBlock(cost_function,
//                                                     p_LossFunction,
//                                                     &map_poses.at(view->id_pose)[0],
//                                                     structure_landmark_it.second.X.data());
                        }
                    }
                    else
                    {
                        std::cerr << "Cannot create a CostFunction for this camera model." << std::endl;
                        return false;
                    }
                }
                if (options.structure_opt == Structure_Parameter_Type::NONE)
                    problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
            }

            if (options.control_point_opt.bUse_control_points)
            {
                std::runtime_error("not define");
            }

            // Configure a BA engine and run it
            //  Make Ceres automatically detect the bundle structure.
            ceres::Solver::Options ceres_config_options;
            ceres_config_options.max_num_iterations = 500;
            ceres_config_options.preconditioner_type =
                    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
            ceres_config_options.linear_solver_type =
                    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
            ceres_config_options.sparse_linear_algebra_library_type =
                    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
//            ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
//            ceres_config_options.logging_type = ceres::SILENT;
            ceres_config_options.num_threads = ceres_options_.nb_threads_;
#if CERES_VERSION_MAJOR < 2
            ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
#endif
            ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

            // Solve BA
            ceres::Solver::Summary summary;
            ceres::Solve(ceres_config_options, &problem, &summary);
            if (ceres_options_.bCeres_summary_)
                std::cout << summary.FullReport() << std::endl;

            // If no error, get back refined parameters
            if (!summary.IsSolutionUsable())
            {
                if (ceres_options_.bVerbose_)
                    std::cout << "Bundle Adjustment failed." << std::endl;
                return false;
            }
            else // Solution is usable
            {
                if (ceres_options_.bVerbose_) {
                    // Display statistics about the minimization
                    std::cout << std::endl
                              << "Bundle Adjustment statistics (approximated RMSE):\n"
                              << " #views: " << sfm_data.views.size() << "\n"
                              << " #poses: " << sfm_data.poses.size() << "\n"
                              << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
                              << " #tracks: " << sfm_data.structure.size() << "\n"
                              << " #residuals: " << summary.num_residuals << "\n"
                              << " Initial RMSE: " << std::sqrt(summary.initial_cost / summary.num_residuals) << "\n"
                              << " Final RMSE: " << std::sqrt(summary.final_cost / summary.num_residuals) << "\n"
                              << " Time (s): " << summary.total_time_in_seconds << "\n"
                              << std::endl;
//                    if (options.use_motion_priors_opt)
//                        std::cout << "Usable motion priors: " << (int) b_usable_prior << std::endl;
                }

                // Update camera poses with refined data
                if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE) {
                    for (auto &pose_it : sfm_data.poses) {
                        const IndexT indexPose = pose_it.first;

                        Eigen::Quaterniond Qwi( map_poses.at(indexPose)[6], map_poses.at(indexPose)[3], map_poses.at(indexPose)[4], map_poses.at(indexPose)[5] );
                        Vec3 twi(map_poses.at(indexPose)[0],
                                 map_poses.at(indexPose)[1],
                                 map_poses.at(indexPose)[2]);
                        Mat3 Rwi = Qwi.toRotationMatrix();

                        Vec3 tci = - sfm_data.IG_Ric.transpose() * sfm_data.IG_tic;
                        Mat3 Rcw = ( Rwi * sfm_data.IG_Ric ).transpose();
                        Vec3 twc = twi + Rwi * sfm_data.IG_tic;
//                        Vec3 tiw = - Rwi.transpose() * twi;
//                        Vec3 tcw =  sfm_data.IG_Ric * tci + tiw;
                        // Update the pose
                        Pose3 &pose = pose_it.second;
                        pose = Pose3(Rcw, twc);
                    }
                }

                // Update camera intrinsics with refined data
                if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE) {
                    for (auto &intrinsic_it : sfm_data.intrinsics) {
                        const IndexT indexCam = intrinsic_it.first;

                        const std::vector<double> &vec_params = map_intrinsics.at(indexCam);
                        intrinsic_it.second->updateFromParams(vec_params);
                    }
                }

                // Structure is already updated directly if needed (no data wrapping)

                return true;
            }
        }

        bool Bundle_Adjustment_IMU_Ceres::Adjust(sfm::SfM_Data &sfm_data, const Optimize_Options &options)
        {
            if( options.use_motion_priors_opt )
            {
                std::cerr << "not define yet" << std::endl;
                assert(0);
            }

            ceres::Problem problem;

            double ex_paparm[7];
            {
                Mat3 Ric = sfm_data.IG_Ric;
                Eigen::Quaterniond Qic(Ric);
                Vec3 tic = sfm_data.IG_tic;
                ex_paparm[0] = tic(0);
                ex_paparm[1] = tic(1);
                ex_paparm[2] = tic(2);
                ex_paparm[3] = Qic.x();
                ex_paparm[4] = Qic.y();
                ex_paparm[5] = Qic.z();
                ex_paparm[6] = Qic.w();

                ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
                problem.AddParameterBlock(ex_paparm, 7, local_parameterization);  // p,q
                problem.SetParameterBlockConstant(ex_paparm);
            }

            // Data wrapper for refinement:
            Hash_Map<IndexT, std::vector<double>> map_intrinsics;
            Hash_Map<IndexT, std::vector<double>> map_poses;
            Hash_Map<IndexT, std::vector<double>> map_speed;

            // Setup Poses data & subparametrization
            for (const auto & pose_it : sfm_data.poses)
            {
                const IndexT indexPose = pose_it.first;

                const Pose3 & pose = pose_it.second;
                const Mat3 Rwc = pose.rotation().transpose();
                const Vec3 twc = pose.center();

                Vec3 tci = - sfm_data.IG_Ric.transpose() * sfm_data.IG_tic;
                Mat3 Rwi = Rwc * sfm_data.IG_Ric.transpose();
                Eigen::Quaterniond Qwi(Rwi);
                Vec3 twi = twc + Rwc * tci;

                // angleAxis + translation
                map_poses[indexPose] = {
                        twi(0),
                        twi(1),
                        twi(2),
                        Qwi.x(),
                        Qwi.y(),
                        Qwi.z(),
                        Qwi.w()
                };

                double * parameter_block = &map_poses.at(indexPose)[0];
                ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
                problem.AddParameterBlock(parameter_block, 7, local_parameterization);  // p,q
//                problem.AddParameterBlock(parameter_block, 6);
                if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
                {
                    // set the whole parameter block as constant for best performance
                    problem.SetParameterBlockConstant(parameter_block);
                }
                else  // Subset parametrization
                {
                    std::vector<int> vec_constant_extrinsic;
                    // If we adjust only the translation, we must set ROTATION as constant
                    if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
                    {
                        // Subset rotation parametrization
                        vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), {0,1,2});
                    }
                    // If we adjust only the rotation, we must set TRANSLATION as constant
                    if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
                    {
                        // Subset translation parametrization
                        vec_constant_extrinsic.insert(vec_constant_extrinsic.end(), {3,4,5,6});
                    }
                    if (!vec_constant_extrinsic.empty())
                    {
                        assert(0);
                        ceres::SubsetParameterization *subset_parameterization =
                                new ceres::SubsetParameterization(7, vec_constant_extrinsic);
//                        auto sssssssson =
//                                new ceres::SubsetParameterization(7, vec_constant_extrinsic, local_parameterization);
                        problem.SetParameterization(parameter_block, subset_parameterization);
                    }
                }
            }

            for( const auto& imu:sfm_data.imus )
            {
                const IndexT indexSpd = imu.first;
                Vec3 speedV3d = sfm_data.Speeds.at(indexSpd).speed_;
                Vec3 Ba = imu.second.linearized_ba_;
                Vec3 Bg = imu.second.linearized_bg_;
                map_speed[indexSpd] = {
                        speedV3d(0),
                        speedV3d(1),
                        speedV3d(2),

                        Ba(0),
                        Ba(1),
                        Ba(2),

                        Bg(0),
                        Bg(1),
                        Bg(2)

                };

                double * parameter_block = &map_speed.at(indexSpd)[0];
                problem.AddParameterBlock(parameter_block, map_speed.at(indexSpd).size());
            }

            // Setup Intrinsics data & subparametrization
            for (const auto & intrinsic_it : sfm_data.intrinsics)
            {
                const IndexT indexCam = intrinsic_it.first;

                if (isValid(intrinsic_it.second->getType()))
                {
                    map_intrinsics[indexCam] = intrinsic_it.second->getParams();
                    if (!map_intrinsics.at(indexCam).empty())
                    {
                        double * parameter_block = &map_intrinsics.at(indexCam)[0];
                        problem.AddParameterBlock(parameter_block, map_intrinsics.at(indexCam).size());
                        if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
                        {
                            // set the whole parameter block as constant for best performance
                            problem.SetParameterBlockConstant(parameter_block);
                        }
                        else
                        {
                            const std::vector<int> vec_constant_intrinsic =
                                    intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
                            if (!vec_constant_intrinsic.empty())
                            {
                                ceres::SubsetParameterization *subset_parameterization =
                                        new ceres::SubsetParameterization(
                                                map_intrinsics.at(indexCam).size(), vec_constant_intrinsic);
                                problem.SetParameterization(parameter_block, subset_parameterization);
                            }
                        }
                    }
                }
                else
                {
                    std::cerr << "Unsupported camera type." << std::endl;
                }
            }

            // Set a LossFunction to be less penalized by false measurements
            //  - set it to nullptr if you don't want use a lossFunction.
            ceres::LossFunction * p_LossFunction =
                    ceres_options_.bUse_loss_function_ ?
                    new ceres::HuberLoss(Square(4.0))
                                                       : nullptr;

            // For all visibility add reprojections errors:
            for (auto & structure_landmark_it : sfm_data.structure)
            {
                const Observations & obs = structure_landmark_it.second.obs;

                for (const auto & obs_it : obs)
                {
                    // Build the residual block corresponding to the track observation:
                    const View * view = sfm_data.views.at(obs_it.first).get();

                    // Each Residual block takes a point and a camera as input and outputs a 2
                    // dimensional residual. Internally, the cost function stores the observed
                    // image location and compares the reprojection against the observation.
//                    ceres::CostFunction* cost_function =
//                            (new ceres::AutoDiffCostFunction
//                                    <ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 6, 3>(
//                                    new ResidualVISUALWithIMUErrorFunctor_Pinhole_Intrinsic_Radial_K3(obs_it.second.x.data())));
//                            assert( sfm_data.intrinsics.at(view->id_intrinsic).get()->getType() == PINHOLE_CAMERA_RADIAL3 );
//                            IntrinsicsToCostFunction(sfm_data.intrinsics.at(view->id_intrinsic).get(),
//                                                     obs_it.second.x);
                    Eigen::Vector2d ob_i = obs_it.second.x;
                    VISfM_Projection* cost_function = new VISfM_Projection( ob_i );

                    if (cost_function)
                    {
                        if (!map_intrinsics.at(view->id_intrinsic).empty())
                        {
                            problem.AddResidualBlock(cost_function,
                                                     p_LossFunction,
                                                     &map_poses.at(view->id_pose)[0],
                                                     ex_paparm,
                                                     &map_intrinsics.at(view->id_intrinsic)[0],
                                                     structure_landmark_it.second.X.data());

//                            problem.AddResidualBlock(cost_function,
//                                                     p_LossFunction,
//                                                     &map_intrinsics.at(view->id_intrinsic)[0],
//                                                     &map_poses.at(view->id_pose)[0],
//                                                     ex_paparm,
//                                                     structure_landmark_it.second.X.data());
                        }
                        else
                        {
                            assert(0);
//                            problem.AddResidualBlock(cost_function,
//                                                     p_LossFunction,
//                                                     &map_poses.at(view->id_pose)[0],
//                                                     structure_landmark_it.second.X.data());
                        }
                    }
                    else
                    {
                        std::cerr << "Cannot create a CostFunction for this camera model." << std::endl;
                        return false;
                    }
                }
                if (options.structure_opt == Structure_Parameter_Type::NONE)
                    problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
            }

            ceres::LossFunction * imu_LossFunction =
                    new ceres::CauchyLoss(Square(2.0));
            {
                // TODO xinli first pose speed
                auto pose_i = sfm_data.poses.begin(); pose_i++;
                auto pose_j = std::next(pose_i);
                for(; pose_j != sfm_data.poses.end(); pose_j++, pose_i++)
                {
                    const IndexT indexPose = pose_j->first;
                    auto imu_ptr = sfm_data.imus.at(indexPose);
                    if( imu_ptr.sum_dt_ > 10.0 ) continue;

                    IMUFactor* imu_factor = new IMUFactor(imu_ptr);
                    problem.AddResidualBlock(imu_factor, imu_LossFunction,
                                             &map_poses.at(pose_i->first)[0],
                                             &map_speed.at(pose_i->first)[0],
                                             &map_poses.at(pose_j->first)[0],
                                             &map_speed.at(pose_j->first)[0]);
                }
            }

            if (options.control_point_opt.bUse_control_points)
            {
                std::runtime_error("not define");
            }

            // Configure a BA engine and run it
            //  Make Ceres automatically detect the bundle structure.
            ceres::Solver::Options ceres_config_options;
            ceres_config_options.max_num_iterations = 500;
            ceres_config_options.preconditioner_type =
                    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
            ceres_config_options.linear_solver_type =
                    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
            ceres_config_options.sparse_linear_algebra_library_type =
                    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
            ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
            ceres_config_options.logging_type = ceres::SILENT;
            ceres_config_options.num_threads = ceres_options_.nb_threads_;
#if CERES_VERSION_MAJOR < 2
            ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
#endif
            ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

            // Solve BA
            ceres::Solver::Summary summary;
            ceres::Solve(ceres_config_options, &problem, &summary);
            if (ceres_options_.bCeres_summary_)
                std::cout << summary.FullReport() << std::endl;

            // If no error, get back refined parameters
            if (!summary.IsSolutionUsable())
            {
                if (ceres_options_.bVerbose_)
                    std::cout << "Bundle Adjustment failed." << std::endl;
                return false;
            }
            else // Solution is usable
            {
                if (ceres_options_.bVerbose_) {
                    // Display statistics about the minimization
                    std::cout << std::endl
                              << "Bundle Adjustment With IMU statistics (approximated RMSE):\n"
                              << " #views: " << sfm_data.views.size() << "\n"
                              << " #poses: " << sfm_data.poses.size() << "\n"
                              << " #points: " << sfm_data.structure.size() << "\n"
                              << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
                              << " #tracks: " << sfm_data.structure.size() << "\n"
                              << " #residuals: " << summary.num_residuals << "\n"
                              << " Initial RMSE: " << std::sqrt(summary.initial_cost / summary.num_residuals) << "\n"
                              << " Final RMSE: " << std::sqrt(summary.final_cost / summary.num_residuals) << "\n"
                              << " Time (s): " << summary.total_time_in_seconds << "\n"
                              << std::endl;
//                    if (options.use_motion_priors_opt)
//                        std::cout << "Usable motion priors: " << (int) b_usable_prior << std::endl;
                }

                // Update camera poses with refined data
                if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE) {
                    for (auto &pose_it : sfm_data.poses) {
                        const IndexT indexPose = pose_it.first;

                        Eigen::Quaterniond Qwi( map_poses.at(indexPose)[6], map_poses.at(indexPose)[3], map_poses.at(indexPose)[4], map_poses.at(indexPose)[5] );
                        Vec3 twi(map_poses.at(indexPose)[0],
                                 map_poses.at(indexPose)[1],
                                 map_poses.at(indexPose)[2]);
                        Mat3 Rwi = Qwi.toRotationMatrix();

                        Vec3 tci = - sfm_data.IG_Ric.transpose() * sfm_data.IG_tic;
                        Mat3 Rcw = ( Rwi * sfm_data.IG_Ric ).transpose();
                        Vec3 twc = twi + Rwi * sfm_data.IG_tic;
//                        Vec3 tiw = - Rwi.transpose() * twi;
//                        Vec3 tcw =  sfm_data.IG_Ric * tci + tiw;
                        // Update the pose
                        pose_it.second.SetRoation(Rcw);
                        pose_it.second.SetCenter(twc);
//                        Pose3 &pose = pose_it.second;
//                        pose = Pose3(Rcw, twc);
                    }
                }

                // Update camera intrinsics with refined data
                if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE) {
                    for (auto &intrinsic_it : sfm_data.intrinsics) {
                        const IndexT indexCam = intrinsic_it.first;

                        const std::vector<double> &vec_params = map_intrinsics.at(indexCam);
                        intrinsic_it.second->updateFromParams(vec_params);
                    }
                }

                for( auto& imu:sfm_data.imus )
                {
                    const IndexT indexIMU = imu.first;
                    Vec3 speed(map_speed.at(indexIMU)[0], map_speed.at(indexIMU)[1], map_speed.at(indexIMU)[2]);
                    Vec3 ba(map_speed.at(indexIMU)[3], map_speed.at(indexIMU)[4], map_speed.at(indexIMU)[5]);
                    Vec3 bg(map_speed.at(indexIMU)[6], map_speed.at(indexIMU)[7], map_speed.at(indexIMU)[8]);

                    sfm_data.Speeds.at(indexIMU).speed_ = speed;
                    imu.second.linearized_ba_ = ba;
                    imu.second.linearized_bg_ = bg;
                }

                // Structure is already updated directly if needed (no data wrapping)

                return true;
            }
        }

    }
}