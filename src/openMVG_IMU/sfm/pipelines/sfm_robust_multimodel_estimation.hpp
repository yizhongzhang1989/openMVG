// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

  struct RelativePose_MultiInfo
{
  Mat3 model_matrix;
  bool b_coplanar;         // false: recover pose from essential , true: recover pose from homography
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
 *  by using a multiple geometry model.
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
  RelativePose_MultiInfo & relativePose_essential,
  RelativePose_MultiInfo & relativePose_homography,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count = 4096,
  int view_id1 = -1,
  int view_id2 = -1
);

/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a robust essential matrix estimation which validated by imu relative pose.
 *
 * @param[in] intrinsics1 camera 1 intrinsics
 * @param[in] intrinsics2 camera 2 intrinsics
 * @param[in] x1 image points in image 1
 * @param[in] x2 image points in image 2
 * @param[out] relativePose_info relative pose information
 * @param[in] relativePose_info relative pose information from imu
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
 *  by using a robust homography matrix estimation.
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
 *  by using a robust essential matrix estimation.
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
