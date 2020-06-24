// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

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

// Decompose an homography matrix into the possible rotations, translations,
// and plane normal vectors, according to:
//
//    Malis, Ezio, and Manuel Vargas. "Deeper understanding of the homography
//    decomposition for vision-based control." (2007): 90.
//
// The first pose is assumed to be P = [I | 0]. Note that the homography is
// plane-induced if `R.size() == t.size() == n.size() == 4`. If `R.size() ==
// t.size() == n.size() == 1` the homography is pure-rotational.
//
// @param[in] H           3x3 homography matrix.
// @param[in] K           3x3 calibration matrix.
// @param[out] R          Possible 3x3 rotation matrices.
// @param[out] t          Possible translation vectors.
// @param[out] n          Possible normal vectors.
// (taken from COLMAP without modification)
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
* (taken from OpenMVG::RelativePoseFromEssential with modification)
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
// (taken from COLMAP with modification)
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
* @param[in] bearing1 Bearing vector of the first camera
* @param[in] bearing2 Bearing vector of the second camera
* @param[in] R rotation of relative pose between first and second camera
* @param[in] t translation of relative pose between first and second camera
* @param[out] points3D_cmb The 3d points in "front" of two cameras.
* (taken from COLMAP with modification)
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
// @param[in] proj_matrix1 3x4 first projection matrix.
// @param[in] proj_matrix2 3x4 second projection matrix.
// @param[in] points3D     input 3d points.
// @param[in] kMinDepth    minimum depth threshold of 3d points.
// @param[in] max_depth    maximum depth thershold of 3d points.
// @return                whether the 3d point lies in front of both cameras.
// (taken from COLMAP with modification)
bool CheiralityTest_COLMAP(const Eigen::Matrix<double,3,4>& proj_matrix1,const Eigen::Matrix<double,3,4>& proj_matrix2,
                           const Eigen::Vector3d& point3D,const double& kMinDepth,const double& max_depth);

// Calculate depth of 3D point with respect to camera.
//
// The depth is defined as the Euclidean distance of a 3D point from the
// camera and is positive if the 3D point is in front and negative if
// behind of the camera.
//
// @param[in] proj_matrix     3x4 projection matrix.
// @param[in] point3D         3D point as 3x1 vector.
//
// @return                Depth of 3D point.
// (taken from COLMAP without modification)
double CalculateDepth(const Eigen::Matrix<double,3,4>& proj_matrix,
                      const Eigen::Vector3d& point3D);

// Compose homography matrix from relative pose.
//
// @param[in] K1      3x3 calibration matrix of first camera.
// @param[in] K2      3x3 calibration matrix of second camera.
// @param[in] R       Most probable 3x3 rotation matrix.
// @param[in] t       Most probable 3x1 translation vector.
// @param[in] n       Most probable 3x1 normal vector.
// @param[in] d       Orthogonal distance from plane.
//
// @return        3x3 homography matrix.
// (taken from COLMAP with modification)
Eigen::Matrix3d HomographyMatrixFromPose(const Eigen::Matrix3d& K1,
                                         const Eigen::Matrix3d& K2,
                                         const Eigen::Matrix3d& R,
                                         const Eigen::Vector3d& t,
                                         const Eigen::Vector3d& n,
                                         const double d);
// select optimal homography matrix from known correspondings.
//
// @param[in] K1       3x3 calibration matrix of first camera.
// @param[in] K2       3x3 calibration matrix of second camera.
// @param[in] points1  First set of corresponding points(image frame).
// @param[in] points2  Second set of corresponding points(image frame).
// @param[in] bearing1 First set of corresponding normalized points(camera frame).
// @param[in] bearing2 Second set of corresponding normalized points(camera frame).
//
// @return        3x3 homography matrix.
// (owned by BC)
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
// @param[in] x1                 Points tested.
// @param[in] n                  3x1 normal vector of plane.
// @param[in] d                  3x1 bias of plane.
// @param[in] maxErrorThreshold  Maximum error thershold.
// @param[out] num_valid          Number of points satisfiying the thershold.
//
// @return        the sum of error of all points.
// (owned by BC)
double CalculateErrorPlane(const Eigen::Matrix<double, 3, Eigen::Dynamic>& x1,
                     const Eigen::Vector3d n,const double d,const double maxErrorThreshold,size_t& num_valid);

// count the number of points satisfying the input homography.
// error = norm(x2 - H * x1)
//
// @param[in] x1                 First set of corresponding points(image frame).
// @param[in] x2                 Second set of corresponding points(image frame).
// @param[in] H                  homography matrix.
// @param[in] maxNormThreshold   Maximum error thershold.
// @param[out] num_valid          Number of points satisfiying the thershold.
//
// @return        the sum of error of all points.
// (owned by BC)
double CalculateErrorH(const Eigen::Matrix<double, 2, Eigen::Dynamic>& x1,
                     const Eigen::Matrix<double, 2, Eigen::Dynamic>& x2,
                     const Eigen::Matrix3d& H,const double maxNormThreshold,size_t& num_valid);

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_HOMOGRAPHY_MATRIX_UTILS_H_
