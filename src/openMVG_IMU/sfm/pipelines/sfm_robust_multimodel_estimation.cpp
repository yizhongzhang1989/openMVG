// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

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
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
)
{
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  //recover pose from essential matrix
  RelativePose_MultiInfo relativePose_essential;
  relativePose_essential.b_coplanar = false;
  
  bool flag_essential = robustRelativePose_Essential(intrinsics1,intrinsics2,
                                x1,x2,relativePose_essential,
                                size_ima1,size_ima2,max_iteration_count);
  
  //recover pose from homography matrix
  RelativePose_MultiInfo relativePose_homography;
  relativePose_homography.b_coplanar = true;
  
  int flag_homography = robustRelativePose_Homography(intrinsics1,intrinsics2,
                                x1,x2,relativePose_homography,
                                size_ima1,size_ima2,max_iteration_count);

  size_t En = relativePose_essential.vec_inliers.size();
  size_t Hn = relativePose_homography.vec_inliers.size();

  //detect degenerate case:
  //  1. if number of inliers of homography matrix is larger than 0.8 times number of inliers of essential matrix,
  //     the initial pair is degenerate and should be discarded.
  //  2. otherwise,if estimation of essential matrix succeeds, the initial pair can recover pose from essential matrix.
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
  //END(Author: BC)===================================================
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
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  relativePose_info.b_coplanar = false;
  //recover pose from essential matrix
  bool flag_essential = robustRelativePose_Essential(intrinsics1,intrinsics2,
                                x1,x2,relativePose_info,
                                size_ima1,size_ima2,max_iteration_count);
  const double MaximumAngularError = 5.0;
  if(flag_essential)
  {
    //calculate the angular error between reltive poses from imu and from essential matrix,
    //and the initial pair is accepted if the error is less than maximum angular error(threshold).
    const Mat3 R1 = imu_relative_pose.rotation(); //imu
    const Mat3 R2 = relativePose_info.relativePose.rotation(); // essential
    const double angularErrorDegree = R2D(getRotationMagnitude(R1 * R2.transpose()));
    if(angularErrorDegree < MaximumAngularError)
    {
      return true;
    }

  }
  return false;
  //END(Author: BC)===================================================
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

    //relativePose_info.found_residual_precision = ac_ransac_output.first;

    if (relativePose_info.vec_inliers.size() <
        2.5 * KernelType::Solver::MINIMUM_SAMPLES )
    {
       //pair failed in inlier threshold
      return 0; // no sufficient coverage (the model does not support enough samples)
    }
  }
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
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
    if(vec_selected_points.size() == 0)
    {
      return 0;  //pair failed in 3d points threshold
    }
    return -1;   //pair failed in ratio threshold 
  }
  relativePose_info.relativePose = relative_pose;
  return 1;  //pair succeed in homography
  //END(Author: BC)===================================================
  

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
