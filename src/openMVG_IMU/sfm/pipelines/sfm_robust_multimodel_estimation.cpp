// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/sfm/pipelines/sfm_robust_multimodel_estimation.hpp"
#include "openMVG_IMU/multiview/homography_matrix.hpp"
#include <array>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include <fstream>

using namespace openMVG::cameras;
using namespace openMVG::geometry;

namespace openMVG {
namespace sfm {

bool robustRelativePose_MultiModel
(
  const IntrinsicBase * intrinsics1,
  const IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  RelativePose_MultiInfo & relativePose_essential,
  RelativePose_MultiInfo & relativePose_homography,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count,
  int view_id1,
  int view_id2
)
{

  
  //RelativePose_MultiInfo relativePose_essential;
  relativePose_essential.b_coplanar = false;
  
  bool flag_essential = robustRelativePose_Essential(intrinsics1,intrinsics2,
                                x1,x2,relativePose_essential,
                                size_ima1,size_ima2,max_iteration_count);
  //RelativePose_MultiInfo relativePose_homography;
  relativePose_homography.b_coplanar = true;
  
  int flag_homography = robustRelativePose_Homography(intrinsics1,intrinsics2,
                                x1,x2,relativePose_homography,
                                size_ima1,size_ima2,max_iteration_count);

  size_t En = relativePose_essential.vec_inliers.size();
  size_t Hn = relativePose_homography.vec_inliers.size();
  
  if((double)Hn > (double)0.8 * En)
  {
    
    return false;
    // if(flag_homography != 0)
    // {
    //     relativePose_info = relativePose_homography;
    // }
    // else
    // {
    //   return false;
    // }
    
  }
  else if(flag_essential)
  {
    relativePose_info = relativePose_essential;
  }
  else
  {
    return false;
  }
    
  
  return true;

}

bool robustRelativePose_IMU(
  const IntrinsicBase * intrinsics1,
  const IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const Pose3& imu_relative_pose,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
  
)
{
  relativePose_info.b_coplanar = false;
  
  bool flag_essential = robustRelativePose_Essential(intrinsics1,intrinsics2,
                                x1,x2,relativePose_info,
                                size_ima1,size_ima2,max_iteration_count);
  const double MaximumAngularError = 5.0;
  if(flag_essential)
  {
    const Mat3 R1 = imu_relative_pose.rotation(); //imu
    const Mat3 R2 = relativePose_info.relativePose.rotation(); // essential
    const double angularErrorDegree = R2D(getRotationMagnitude(R1 * R2.transpose()));
    if(angularErrorDegree < MaximumAngularError)
    {
      return true;
    }

  }
  return false;
  
}


int robustRelativePose_Homography
(
  const IntrinsicBase * intrinsics1,
  const IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
)
{
  if (!intrinsics1 || !intrinsics2)
    return false;

  // Compute the bearing vectors
  const Mat3X
    bearing1 = (*intrinsics1)(x1),
    bearing2 = (*intrinsics2)(x2);

  
  {
    // Define the AContrario adaptor to use the 5 point essential matrix solver.
    using KernelType =
      openMVG::robust::ACKernelAdaptor<
        openMVG::homography::kernel::FourPointSolver,
        openMVG::homography::kernel::AsymmetricError,
        UnnormalizerI,
        Mat3>;
    KernelType kernel(
      x1, size_ima1.first, size_ima1.second,
      x2, size_ima2.first, size_ima2.second,
      false); // configure as point to point error model.
    // The Model type

    // Call the Robust Estimator on the KERNEL
    const std::pair<double,double> ACRansacOut = // Return the precision & the associated NFA
      ACRANSAC(
        kernel, // The Kernel (embed the correspondences, the Model Solver & the fitting metric.
        relativePose_info.vec_inliers, // The inlier list
        max_iteration_count, // Max iteration count (it can be less if a meaningful is found at the first iterations)
        &relativePose_info.model_matrix, // Returned model
        relativePose_info.initial_residual_tolerance, // No apriori Threshold
        false // Verbose to console
       );
    relativePose_info.found_residual_precision = ACRansacOut.first;

    ////////bc debug start /////////
    bool bdebug = false;
    if(bdebug)
    {
      std::ofstream homography_matrix_log("D:/Microsoft_internship/data/synthetic_data/pharmacy/key_openMVG_seq_recon_of_filtering_phase/homography_error.txt");
      const Mat3 H = relativePose_info.model_matrix;
      std::vector<double> norm_s;
      for(int i = 0;i<x1.cols();i++)
      {
        homography_matrix_log<<i+1<<":\n";
        Vec3 p1(x1.col(i)(0),x1.col(i)(1),1.0);
        Vec3 p2(x2.col(i)(0),x2.col(i)(1),1.0);
        homography_matrix_log <<"p1:"<<p1(0)<<" "<<p1(1)<<" "<<p1(2)<<"\n";
        homography_matrix_log <<"p2:"<<p2(0)<<" "<<p2(1)<<" "<<p2(2)<<"\n";
        Vec3 Hp1 = H*p1;
        Hp1 /= Hp1(2);
        homography_matrix_log <<"Hp1:"<<Hp1(0)<<" "<<Hp1(1)<<" "<<Hp1(2)<<"\n";
        double norm_ = (Hp1 - p2).norm();
        norm_s.push_back(norm_);
        homography_matrix_log <<"norm:"<<norm_<<"\n";
      }
      std::sort(norm_s.begin(),norm_s.end());
      homography_matrix_log <<"norm_s:\n";
      for(int i = 0 ;i < norm_s.size();i++)
      {
        homography_matrix_log<<norm_s[i]<<"\n";
      }
      homography_matrix_log.close();
    }
    ////////bc debug end /////////

    //relativePose_info.found_residual_precision = ac_ransac_output.first;

    if (relativePose_info.vec_inliers.size() <
        2.5 * KernelType::Solver::MINIMUM_SAMPLES )
    {
       //pair failed in inlier threshold, which should be discard
      return 0; // no sufficient coverage (the model does not support enough samples)
    }
  }
  // Eigen::Matrix3d R;
  // Eigen::Vector3d t,n;
  // std::vector<Eigen::Vector3d> points3D;
  // bool flag = PoseFromHomographyMatrix(relativePose_info.model_matrix,
  //                     dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics1)->K(),
  //                     dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics2)->K(),
  //                             //const std::vector<Eigen::Vector2d>& points1,
  //                             //const std::vector<Eigen::Vector2d>& points2,
  //                             x1,
  //                             x2,
  //                             bearing1,
  //                             bearing2,
  //                             &R, &t,
  //                             &n,
  //                             &points3D);
  // relativePose_info.relativePose = Pose3(R , -R.transpose() * t);
  // return flag;
  
   // estimation of the relative poses based on the cheirality test
  Pose3 relative_pose;
  std::vector<uint32_t> vec_selected_points;
  if (!RelativePoseFromHomography(
    relativePose_info.model_matrix,
    dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics1)->K(),
    dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics2)->K(),
    bearing1,
    bearing2,

    relativePose_info.vec_inliers, &relative_pose,
    &vec_selected_points))
  {
    //relativePose_info.relativePose = relative_pose;   //////6/16 must deleted
    if(vec_selected_points.size() == 0)
    {
      return 0;  //pair failed in 3d points threshold, which should be discard
    }
    return -1;   //pair failed in ratio threshold , which should be discard
  }
  relativePose_info.relativePose = relative_pose;
  return 1;  //pair succeed in homography , which can recover from essential or from homography

  

}

bool robustRelativePose_Essential
(
  const IntrinsicBase * intrinsics1,
  const IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
)
{
  if (!intrinsics1 || !intrinsics2)
    return false;

  // Compute the bearing vectors
  const Mat3X
    bearing1 = (*intrinsics1)(x1),
    bearing2 = (*intrinsics2)(x2);

  if (isPinhole(intrinsics1->getType())
      && isPinhole(intrinsics2->getType()))
  {
    // Define the AContrario adaptor to use the 5 point essential matrix solver.
    using KernelType = robust::ACKernelAdaptorEssential<
      openMVG::essential::kernel::FivePointSolver,
      openMVG::fundamental::kernel::EpipolarDistanceError,
      Mat3>;
    KernelType kernel(x1, bearing1, size_ima1.first, size_ima1.second,
                      x2, bearing2, size_ima2.first, size_ima2.second,
                      dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics1)->K(),
                      dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics2)->K());

    // Robustly estimation of the Model and its precision
    const auto ac_ransac_output = robust::ACRANSAC(
      kernel, relativePose_info.vec_inliers,
      max_iteration_count, &relativePose_info.model_matrix,
      relativePose_info.initial_residual_tolerance, false);

    relativePose_info.found_residual_precision = ac_ransac_output.first;

    if (relativePose_info.vec_inliers.size() <
        2.5 * KernelType::Solver::MINIMUM_SAMPLES )
    {
      return false; // no sufficient coverage (the model does not support enough samples)
    }
  }
  else
  {
    // Define the AContrario adaptor to use the 8 point essential matrix solver.
    typedef openMVG::robust::ACKernelAdaptor_AngularRadianError<
        openMVG::EightPointRelativePoseSolver,
        // openMVG::essential::kernel::FivePointSolver,
        openMVG::AngularError,
        Mat3>
        KernelType;

    KernelType kernel(bearing1, bearing2);

    // Robustly estimate the Essential matrix with A Contrario ransac
    const double upper_bound_precision =
      (relativePose_info.initial_residual_tolerance == std::numeric_limits<double>::infinity()) ?
        std::numeric_limits<double>::infinity()
        : D2R(relativePose_info.initial_residual_tolerance);
    const auto ac_ransac_output =
      ACRANSAC(kernel, relativePose_info.vec_inliers,
        max_iteration_count, &relativePose_info.model_matrix,
        upper_bound_precision, false);

    const double & threshold = ac_ransac_output.first;
    relativePose_info.found_residual_precision = R2D(threshold); // Degree

    if (relativePose_info.vec_inliers.size() <
        2.5 * KernelType::Solver::MINIMUM_SAMPLES )
    {
      return false; // no sufficient coverage (the model does not support enough samples)
    }
  }

  // estimation of the relative poses based on the cheirality test
  Pose3 relative_pose;
  if (!RelativePoseFromEssential(
    bearing1,
    bearing2,
    relativePose_info.model_matrix,
    relativePose_info.vec_inliers, &relative_pose))
  {
    return false;
  }
  relativePose_info.relativePose = relative_pose;
  return true;
}

} // namespace sfm
} // namespace openMVG
