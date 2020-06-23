// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#ifndef COLMAP_SRC_BASE_HOMOGRAPHY_MATRIX_UTILS_H_
#define COLMAP_SRC_BASE_HOMOGRAPHY_MATRIX_UTILS_H_

#include <vector>

#include <Eigen/Core>

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/multiview/triangulation_method.hpp"


namespace openMVG {

template <typename T>
int SignOfNumber(const T val) {
  return (T(0) < val) - (val < T(0));
}

// (COLMAP)Decompose an homography matrix into the possible rotations, translations,
// and plane normal vectors, according to:
//
//    Malis, Ezio, and Manuel Vargas. "Deeper understanding of the homography
//    decomposition for vision-based control." (2007): 90.
//
// The first pose is assumed to be P = [I | 0]. Note that the homography is
// plane-induced if `R.size() == t.size() == n.size() == 4`. If `R.size() ==
// t.size() == n.size() == 1` the homography is pure-rotational.
//
// @param H          3x3 homography matrix.
// @param K          3x3 calibration matrix.
// @param R          Possible 3x3 rotation matrices.
// @param t          Possible translation vectors.
// @param n          Possible normal vectors.
void DecomposeHomographyMatrix(const Eigen::Matrix3d& H,
                               const Eigen::Matrix3d& K1,
                               const Eigen::Matrix3d& K2,
                               std::vector<Eigen::Matrix3d>* R,
                               std::vector<Eigen::Vector3d>* t,
                               std::vector<Eigen::Vector3d>* n);
                        
/**
* @brief Estimate the best possible relative pose from H imitating recovery from E in OpenMVG.
*  Four relative poses can be build from the Hmatrix decomposition.
*  We keep the one with most of the point in front of the camera.
*
* @param[in] H homography matrix
* @param[in] K1 3x3 calibration matrix of first camera.
* @param[in] K2 3x3 calibration matrix of second camera.
* @param[in] x1 bearing vectors corresponding to image observation in image 1
* @param[in] x2 bearing vectors corresponding to image observation in image 2
* @param[in] bearing_vector_index_to_use selection of x1, x2 columns that are used
* @param[out] relative_pose the estimated relative pose
* @param[out] vec_selected_points return the index of bearing_vector_index_to_use that are
*    in front of the cameras
* @param[out] vec_points return the 3D point corresponding to
*    vec_selected_points indexes
* @param[in] positive_depth_solution_ratio Pourcentage ratio threshold used
*    to discard if there is two good solution that have many points in front
*    of the cameras
*/

bool RelativePoseFromHomography
(
  
  const Mat3 & H,
  const Eigen::Matrix3d& K1,
  const Eigen::Matrix3d& K2,
  const Mat3X & x1,
  const Mat3X & x2,
  const std::vector<uint32_t> & bearing_vector_index_to_use,
  geometry::Pose3 * relative_pose = nullptr,
  std::vector<uint32_t> * vec_selected_points = nullptr,
  std::vector<Vec3> * vec_points = nullptr,
  const double positive_depth_solution_ratio = 0.7,
  const ETriangulationMethod triangulation_method = ETriangulationMethod::DEFAULT
);

// Recover the most probable pose from the given homography matrix recovery from h in COLMAP.
//
// The pose of the first image is assumed to be P = [I | 0].
//
// @param[in] H            3x3 homography matrix.
// @param[in] K1           3x3 calibration matrix of first camera.
// @param[in] K2           3x3 calibration matrix of second camera.
// @param[in] bearing1     First set of corresponding normalized points(camera frame).
// @param[in] bearing2     Second set of corresponding normalized points(camera frame).

// @param[out] R           Most probable 3x3 rotation matrix.
// @param[out] t           Most probable 3x1 translation vector.
// @param[out] n           Most probable 3x1 normal vector.
// @param[out] points3D    Triangulated 3D points infront of camera
//                         (only if homography is not pure-rotational).
bool PoseFromHomographyMatrix(const Eigen::Matrix3d& H,
                              const Eigen::Matrix3d& K1,
                              const Eigen::Matrix3d& K2,

                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing1,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing2,
                              Eigen::Matrix3d* R, Eigen::Vector3d* t,
                              Eigen::Vector3d* n,
                              std::vector<Eigen::Vector3d>* points3D);

/**
* @brief Test if a 3d point is in "front" of two camera.
*        This function test if a 3d point can be visible according two camera pose
*        and two bearing vector taken from COLMAP.
*
* @param bearing1 Bearing vector of the first camera
* @param bearing2 Bearing vector of the second camera
* @param R rotation of relative pose between first and second camera
* @param t translation of relative pose between first and second camera
* @param points3D_cmb The 3d points in "front" of two cameras.
*
*/
void CheckCheirality(const Eigen::Matrix3d& R,const Eigen::Vector3d& t,
                    const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing1,
                    const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing2,
                    std::vector<Eigen::Vector3d>& points3D_cmb
                    );

// Perform cheirality constraint test, i.e., determine whether the 3d point
// lies in front of both cameras.
// 
//
// @param proj_matrix1 3x4 first projection matrix.
// @param proj_matrix2 3x4 second projection matrix.
// @param points3D     input 3d points.
// @param kMinDepth    minimum depth threshold of 3d points.
// @param max_depth    maximum depth thershold of 3d points.
bool CheiralityTest_COLMAP(const Eigen::Matrix<double,3,4>& proj_matrix1,const Eigen::Matrix<double,3,4>& proj_matrix2,
                           const Eigen::Vector3d& point3D,const double& kMinDepth,const double& max_depth);

// Calculate depth of 3D point with respect to camera.
//
// The depth is defined as the Euclidean distance of a 3D point from the
// camera and is positive if the 3D point is in front and negative if
// behind of the camera.
//
// @param proj_matrix     3x4 projection matrix.
// @param point3D         3D point as 3x1 vector.
//
// @return                Depth of 3D point.
double CalculateDepth(const Eigen::Matrix<double,3,4>& proj_matrix,
                      const Eigen::Vector3d& point3D);

// Compose homography matrix from relative pose.
//
// @param K1      3x3 calibration matrix of first camera.
// @param K2      3x3 calibration matrix of second camera.
// @param R       Most probable 3x3 rotation matrix.
// @param t       Most probable 3x1 translation vector.
// @param n       Most probable 3x1 normal vector.
// @param d       Orthogonal distance from plane.
//
// @return        3x3 homography matrix.
Eigen::Matrix3d HomographyMatrixFromPose(const Eigen::Matrix3d& K1,
                                         const Eigen::Matrix3d& K2,
                                         const Eigen::Matrix3d& R,
                                         const Eigen::Vector3d& t,
                                         const Eigen::Vector3d& n,
                                         const double d);
// select optimal homography matrix from known correspondings.
//
// @param K1       3x3 calibration matrix of first camera.
// @param K2       3x3 calibration matrix of second camera.
// @param points1  First set of corresponding points(image frame).
// @param points2  Second set of corresponding points(image frame).
// @param bearing1 First set of corresponding normalized points(camera frame).
// @param bearing2 Second set of corresponding normalized points(camera frame).
//
// @return        3x3 homography matrix.
Eigen::Matrix3d ComputeHomographyMatrixfromKnownPose(
                              const Eigen::Matrix3d& K1,
                              const Eigen::Matrix3d& K2,
                              const Eigen::Matrix<double, 2, Eigen::Dynamic>& points1,
                              const Eigen::Matrix<double, 2, Eigen::Dynamic>& points2,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing1,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing2,
                              const Eigen::Matrix3d& R,
                              const Eigen::Vector3d& t
                    );
// count the number of points satisfying the input plane.
// error = abs(n.dot(x1) + d)
//
// @param x1                 Points tested.
// @param n                  3x1 normal vector of plane.
// @param d                  3x1 bias of plane.
// @param maxErrorThreshold  Maximum error thershold.
// @param num_valid          Number of points satisfiying the thershold.
//
// @return        the sum of error of all points.
double CalculateErrorPlane(const Eigen::Matrix<double, 3, Eigen::Dynamic>& x1,
                     const Eigen::Vector3d n,const double d,const double maxErrorThreshold,size_t& num_valid);

// count the number of points satisfying the input homography.
// error = norm(x2 - H * x1)
//
// @param x1                 First set of corresponding points(image frame).
// @param x2                 Second set of corresponding points(image frame).
// @param H                  homography matrix.
// @param maxNormThreshold   Maximum error thershold.
// @param num_valid          Number of points satisfiying the thershold.
//
// @return        the sum of error of all points.
double CalculateErrorH(const Eigen::Matrix<double, 2, Eigen::Dynamic>& x1,
                     const Eigen::Matrix<double, 2, Eigen::Dynamic>& x2,
                     const Eigen::Matrix3d& H,const double maxNormThreshold,size_t& num_valid);

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_HOMOGRAPHY_MATRIX_UTILS_H_
