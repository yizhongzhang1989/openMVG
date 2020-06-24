// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_SFM_SFM_ROBUST_MULTIMODEL_ESTIMATION_HPP
#define OPENMVG_IMU_SFM_SFM_ROBUST_MULTIMODEL_ESTIMATION_HPP

#include <limits>
#include <utility>
#include <vector>

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG { namespace cameras { struct IntrinsicBase; } }

namespace openMVG {
namespace sfm {

//(taken from OpenMVG with modification)
  struct RelativePose_MultiInfo
{
  Mat3 model_matrix;
  bool b_coplanar;         //true: the initial pair is a degenerate case
  geometry::Pose3 relativePose;
  std::vector<uint32_t> vec_inliers;
  double initial_residual_tolerance;
  double found_residual_precision;

  RelativePose_MultiInfo()
    :initial_residual_tolerance(std::numeric_limits<double>::infinity()),
    found_residual_precision(std::numeric_limits<double>::infinity()),
    b_coplanar(false)
  {}
};

/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a multiple geometry model(owned by BC).
 *
 * @param[in] intrinsics1 camera 1 intrinsics
 * @param[in] intrinsics2 camera 2 intrinsics
 * @param[in] x1 image points in image 1
 * @param[in] x2 image points in image 2
 * @param[out] relativePose_info relative pose information
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[in] max iteration count
 */
bool robustRelativePose_MultiModel
(
  const cameras::IntrinsicBase * intrinsics1,
  const cameras::IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count = 4096
);

/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a robust essential matrix estimation with imu validation(owned by BC).
 *
 * @param[in] intrinsics1 camera 1 intrinsics
 * @param[in] intrinsics2 camera 2 intrinsics
 * @param[in] x1 image points in image 1
 * @param[in] x2 image points in image 2
 * @param[out] relativePose_info relative pose information
 * @param[in] imu_relative_pose relative pose information from imu
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[in] max iteration count
 */
bool robustRelativePose_IMU(
  const cameras::IntrinsicBase * intrinsics1,
  const cameras::IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const geometry::Pose3& imu_relative_pose,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
  
);
/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a robust homography matrix estimation(taken from OpenMVG with modification).
 *
 * @param[in] intrinsics1 camera 1 intrinsics
 * @param[in] intrinsics2 camera 2 intrinsics
 * @param[in] x1 image points in image 1
 * @param[in] x2 image points in image 2
 * @param[out] relativePose_info relative pose information
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[in] max iteration count
 */
int robustRelativePose_Homography
(
  const cameras::IntrinsicBase * intrinsics1,
  const cameras::IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count = 4096
);
/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a robust essential matrix estimation(taken from OpenMVG::robustRelativePose() without modification).
 *
 * @param[in] intrinsics1 camera 1 intrinsics
 * @param[in] intrinsics2 camera 2 intrinsics
 * @param[in] x1 image points in image 1
 * @param[in] x2 image points in image 2
 * @param[out] relativePose_info relative pose information
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[in] max iteration count
 */
bool robustRelativePose_Essential
(
  const cameras::IntrinsicBase * intrinsics1,
  const cameras::IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_MultiInfo & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count = 4096
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_HPP
