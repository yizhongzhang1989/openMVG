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


#include "openMVG_IMU/multiview/homography_matrix.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include <array>
#include <fstream>  //bc debug
#include <Eigen/Dense>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


namespace  openMVG{
namespace {

double ComputeOppositeOfMinor(const Eigen::Matrix3d& matrix, const size_t row,
                              const size_t col) {
  const size_t col1 = col == 0 ? 1 : 0;
  const size_t col2 = col == 2 ? 1 : 2;
  const size_t row1 = row == 0 ? 1 : 0;
  const size_t row2 = row == 2 ? 1 : 2;
  return (matrix(row1, col2) * matrix(row2, col1) -
          matrix(row1, col1) * matrix(row2, col2));
}

Eigen::Matrix3d ComputeHomographyRotation(const Eigen::Matrix3d& H_normalized,
                                          const Eigen::Vector3d& tstar,
                                          const Eigen::Vector3d& n,
                                          const double v) {
  return H_normalized *
         (Eigen::Matrix3d::Identity() - (2.0 / v) * tstar * n.transpose());
}

}  // namespace



void DecomposeHomographyMatrix(const Eigen::Matrix3d& H,
                               const Eigen::Matrix3d& K1,
                               const Eigen::Matrix3d& K2,
                               std::vector<Eigen::Matrix3d>* R,
                               std::vector<Eigen::Vector3d>* t,
                               std::vector<Eigen::Vector3d>* n) {
  // Remove calibration from homography.
  Eigen::Matrix3d H_normalized = K2.inverse() * H * K1;

  // Remove scale from normalized homography.
  Eigen::JacobiSVD<Eigen::Matrix3d> hmatrix_norm_svd(H_normalized);
  H_normalized.array() /= hmatrix_norm_svd.singularValues()[1];

  const Eigen::Matrix3d S =
      H_normalized.transpose() * H_normalized - Eigen::Matrix3d::Identity();

  // Check if H is rotation matrix.
  const double kMinInfinityNorm = 1e-3;
  if (S.lpNorm<Eigen::Infinity>() < kMinInfinityNorm) {
    *R = {H_normalized};
    *t = {Eigen::Vector3d::Zero()};
    *n = {Eigen::Vector3d::Zero()};
    return;
  }

  const double M00 = ComputeOppositeOfMinor(S, 0, 0);
  const double M11 = ComputeOppositeOfMinor(S, 1, 1);
  const double M22 = ComputeOppositeOfMinor(S, 2, 2);

  const double rtM00 = std::sqrt(M00);
  const double rtM11 = std::sqrt(M11);
  const double rtM22 = std::sqrt(M22);

  const double M01 = ComputeOppositeOfMinor(S, 0, 1);
  const double M12 = ComputeOppositeOfMinor(S, 1, 2);
  const double M02 = ComputeOppositeOfMinor(S, 0, 2);

  const int e12 = SignOfNumber(M12);
  const int e02 = SignOfNumber(M02);
  const int e01 = SignOfNumber(M01);

  const double nS00 = std::abs(S(0, 0));
  const double nS11 = std::abs(S(1, 1));
  const double nS22 = std::abs(S(2, 2));

  const std::array<double, 3> nS{{nS00, nS11, nS22}};
  const size_t idx =
      std::distance(nS.begin(), std::max_element(nS.begin(), nS.end()));

  Eigen::Vector3d np1;
  Eigen::Vector3d np2;
  if (idx == 0) {
    np1[0] = S(0, 0);
    np2[0] = S(0, 0);
    np1[1] = S(0, 1) + rtM22;
    np2[1] = S(0, 1) - rtM22;
    np1[2] = S(0, 2) + e12 * rtM11;
    np2[2] = S(0, 2) - e12 * rtM11;
  } else if (idx == 1) {
    np1[0] = S(0, 1) + rtM22;
    np2[0] = S(0, 1) - rtM22;
    np1[1] = S(1, 1);
    np2[1] = S(1, 1);
    np1[2] = S(1, 2) - e02 * rtM00;
    np2[2] = S(1, 2) + e02 * rtM00;
  } else if (idx == 2) {
    np1[0] = S(0, 2) + e01 * rtM11;
    np2[0] = S(0, 2) - e01 * rtM11;
    np1[1] = S(1, 2) + rtM00;
    np2[1] = S(1, 2) - rtM00;
    np1[2] = S(2, 2);
    np2[2] = S(2, 2);
  }

  const double traceS = S.trace();
  const double v = 2.0 * std::sqrt(1.0 + traceS - M00 - M11 - M22);

  const double ESii = SignOfNumber(S(idx, idx));
  const double r_2 = 2 + traceS + v;
  const double nt_2 = 2 + traceS - v;

  const double r = std::sqrt(r_2);
  const double n_t = std::sqrt(nt_2);

  const Eigen::Vector3d n1 = np1.normalized();
  const Eigen::Vector3d n2 = np2.normalized();

  const double half_nt = 0.5 * n_t;
  const double esii_t_r = ESii * r;

  const Eigen::Vector3d t1_star = half_nt * (esii_t_r * n2 - n_t * n1);
  const Eigen::Vector3d t2_star = half_nt * (esii_t_r * n1 - n_t * n2);

  const Eigen::Matrix3d R1 =
      ComputeHomographyRotation(H_normalized, t1_star, n1, v);
  const Eigen::Vector3d t1 = R1 * t1_star;

  const Eigen::Matrix3d R2 =
      ComputeHomographyRotation(H_normalized, t2_star, n2, v);
  const Eigen::Vector3d t2 = R2 * t2_star;

  *R = {R1, R1, R2, R2};
  *t = {t1, -t1, t2, -t2};
  *n = {-n1, n1, -n2, n2};
}

bool RelativePoseFromHomography
(
  
  const Mat3 & H,
  const Eigen::Matrix3d& K1,
  const Eigen::Matrix3d& K2,
  const Mat3X & x1,
  const Mat3X & x2,
  
  const std::vector<uint32_t> & bearing_vector_index_to_use,
  geometry::Pose3 * relative_pose,
  std::vector<uint32_t> * vec_selected_points,
  std::vector<Vec3> * vec_points,
  const double positive_depth_solution_ratio,
  const ETriangulationMethod triangulation_method
)
{
  using namespace geometry;
  // Recover plausible relative poses from H.
  std::vector<geometry::Pose3> relative_poses;
  assert(x1.cols() == x2.cols());

  std::vector<Eigen::Matrix3d> R_cmbs;
  std::vector<Eigen::Vector3d> t_cmbs;
  std::vector<Eigen::Vector3d> n_cmbs;
  DecomposeHomographyMatrix(H, K1, K2, &R_cmbs, &t_cmbs, &n_cmbs);
  //package the pose
  for(size_t i = 0 ; i < R_cmbs.size() ; i++)
  {
    relative_poses.push_back(geometry::Pose3(R_cmbs[i], -R_cmbs[i].transpose() * t_cmbs[i]));
  }
  
  // Accumulator to find the best solution
  std::vector<uint32_t> cheirality_accumulator(relative_poses.size(), 0);

  // Find which solution is the best:
  // - count how many triangulated observations are in front of the cameras
  std::vector<std::vector<uint32_t>> vec_newInliers(relative_poses.size());
  std::vector<std::vector<Vec3>> vec_3D(relative_poses.size());

  const Pose3 pose1(Mat3::Identity(), Vec3::Zero());

  for (size_t i = 0; i < relative_poses.size(); ++i)
  {
    const Pose3 &pose2 = relative_poses[i];
    Vec3 X;
    for (const uint32_t inlier_idx : bearing_vector_index_to_use)
    {
      const auto
        f1 = x1.col(inlier_idx),
        f2 = x2.col(inlier_idx);
      if (Triangulate2View
      (
        pose1.rotation(), pose1.translation(), f1,
        pose2.rotation(), pose2.translation(), f2,
        X,
        triangulation_method
      ))
      {
        ++cheirality_accumulator[i];
        vec_newInliers[i].push_back(inlier_idx);
        vec_3D[i].push_back(X);
      }
    }
    
  }

  // Check if there is a valid solution:
  const auto iter = std::max_element(cheirality_accumulator.cbegin(), cheirality_accumulator.cend());
  if (*iter == 0)
  {
    // There is no right solution with points in front of the cameras
    return false;
  }

  // Export the best solution data
  const size_t index = std::distance(cheirality_accumulator.cbegin(), iter);
  if (relative_pose)
  {
    (*relative_pose) = relative_poses[index];
  }
  if (vec_selected_points)
    (*vec_selected_points) = vec_newInliers[index];
  if (vec_points)
    (*vec_points) = vec_3D[index];

  // Test if the best solution is good by using the ratio of the two best solution score
  std::sort(cheirality_accumulator.begin(), cheirality_accumulator.end());
  const double ratio = cheirality_accumulator.rbegin()[1]
    / static_cast<double>(cheirality_accumulator.rbegin()[0]);
  return (ratio < positive_depth_solution_ratio);
}


bool PoseFromHomographyMatrix(const Eigen::Matrix3d& H,
                              const Eigen::Matrix3d& K1,
                              const Eigen::Matrix3d& K2,

                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing1,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing2,
                              Eigen::Matrix3d* R, Eigen::Vector3d* t,
                              Eigen::Vector3d* n,
                              std::vector<Eigen::Vector3d>* points3D) {
  
  assert(bearing1.cols() == bearing2.cols());

  // Recover plausible relative poses from H.
  std::vector<Eigen::Matrix3d> R_cmbs;
  std::vector<Eigen::Vector3d> t_cmbs;
  std::vector<Eigen::Vector3d> n_cmbs;
  DecomposeHomographyMatrix(H, K1, K2, &R_cmbs, &t_cmbs, &n_cmbs);
  points3D->clear();

  // Find which solution is the best:
  // - count how many triangulated observations are in front of the cameras
  int num_Maximum_pose = 0;
  for (size_t i = 0; i < R_cmbs.size(); ++i) {
    std::vector<Eigen::Vector3d> points3D_cmb;

    CheckCheirality(R_cmbs[i], t_cmbs[i], bearing1, bearing2, points3D_cmb);
    if (points3D_cmb.size() > points3D->size()) {
      *R = R_cmbs[i];
      *t = t_cmbs[i];
      *n = n_cmbs[i];
      *points3D = points3D_cmb;
      num_Maximum_pose = 1;
    }
    else if(points3D_cmb.size() == points3D->size())
    {
      num_Maximum_pose ++;
    }
  }
  // discard the initial pair that can not be distinguished(i.e. there are multiple best solutions).
  return num_Maximum_pose <= 1;
  //return true;
  
}

void CheckCheirality(const Eigen::Matrix3d& R,const Eigen::Vector3d& t,
                    const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing1,
                    const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing2,
                    std::vector<Eigen::Vector3d>& points3D_cmb
                    )       
{
    assert(bearing1.cols() == bearing2.cols());
    const geometry::Pose3 & pose_1 = geometry::Pose3(Mat3::Identity(), Vec3::Zero());
    const geometry::Pose3 & pose_2 = geometry::Pose3(R,-R.transpose() * t);   //Pose3(R,c)
    //const double kMinDepth = std::numeric_limits<double>::epsilon();
    //const double max_depth = 1000.0f * (R.transpose() * t).norm();
    for(size_t i  = 0 ; i < bearing1.cols() ; i++)
    {
        Vec3 X;
        if(Triangulate2View(pose_1.rotation(),pose_1.translation(),bearing1.col(i),
        pose_2.rotation(), pose_2.translation() ,bearing2.col(i), X)
        )
        {
          if(cameras::CheiralityTest(bearing1.col(i),pose_1,bearing2.col(i),pose_2,X))
          {
              points3D_cmb.push_back(X);
          }

          // if(CheiralityTest_COLMAP(pose_1.asMatrix(),pose_2.asMatrix(),X,kMinDepth,max_depth))
          // {
          //     points3D_cmb_colmap.push_back(X);
          // }
            
        }
    }
}



bool CheiralityTest_COLMAP(const Eigen::Matrix<double,3,4>& proj_matrix1,const Eigen::Matrix<double,3,4>& proj_matrix2,
                           const Eigen::Vector3d& point3D,const double& kMinDepth,const double& max_depth)
{
    const double depth1 = CalculateDepth(proj_matrix1, point3D);
    if (depth1 > kMinDepth && depth1 < max_depth) {
      const double depth2 = CalculateDepth(proj_matrix2, point3D);
      if (depth2 > kMinDepth && depth2 < max_depth) {
        return true;
      }
    }
    return false;
}


double CalculateDepth(const Eigen::Matrix<double,3,4>& proj_matrix,
                      const Eigen::Vector3d& point3D) {
  const double proj_z = proj_matrix.row(2).dot(point3D.homogeneous());
  return proj_z * proj_matrix.col(2).norm();
}

Eigen::Matrix3d HomographyMatrixFromPose(const Eigen::Matrix3d& K1,
                                         const Eigen::Matrix3d& K2,
                                         const Eigen::Matrix3d& R,
                                         const Eigen::Vector3d& t,
                                         const Eigen::Vector3d& n,
                                         const double d) {
  assert(d>=0);
  return K2 * (R - t * n.transpose() / d) * K1.inverse();
}

Eigen::Matrix3d ComputeHomographyMatrixfromKnownPose(
                              const Eigen::Matrix3d& K1,
                              const Eigen::Matrix3d& K2,
                              const Eigen::Matrix<double, 2, Eigen::Dynamic>& points1,
                              const Eigen::Matrix<double, 2, Eigen::Dynamic>& points2,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing1,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& bearing2,
                              const Eigen::Matrix3d& R,
                              const Eigen::Vector3d& t
                    )
{
  ///bearing : the normalized coordinate of features(camera frame)
  ///points: the 2d coordinate of features(image)
    /* initialize random seed: */
    //srand (time(NULL));
    
    //////////////////////////////////////////////////
    ///////////calculate the best H ////////////////
    //////////////////////////////////////////////////
    Vec3 a,b,c;
    double minNumNorm = std::numeric_limits<double>::max();
    Eigen::Matrix3d bestH;
    //find which Hmatrix is best in all plausible Hmatrix:
    //   -calculate the sum of errors according formula(x2 = H * x1)
    for(size_t i = 0; i < bearing1.cols(); i++)
    {
      a = bearing1.col(i);
      for(size_t j = i + 1 ; j < bearing1.cols() ; j++)
      {
        b = bearing1.col(j);
        for(size_t k = j + 1; k< bearing1.cols() ; k++)
        {

        c = bearing1.col(k);

        Vec3 norm = (a-b).cross((c-b));
        //norm = norm.normalized();
        double d= -(norm.dot(a));
        if(d < 0)
        {
          norm = -norm;
          d = -d;
        }
        Eigen::Matrix3d H = HomographyMatrixFromPose(K1,
                                    K2,
                                    R,
                                    t,
                                    norm,
                                    d);
        size_t num_validnorm;
        double norm_ = CalculateErrorH(points1,points2,H,100.0,num_validnorm);
        size_t num_validerror;
        double error_ = CalculateErrorPlane(bearing1,norm,d,0.002,num_validerror);
        
        if(norm_ < minNumNorm)
        {
          minNumNorm = norm_;
          bestH = H;
        }
        }
      }
    }
    
    return bestH;


}

double CalculateErrorPlane(const Eigen::Matrix<double, 3, Eigen::Dynamic>& x1,
                     const Eigen::Vector3d n,const double d,const double maxErrorThreshold,size_t& num_valid)
{
  num_valid = 0;
  double error_sum = 0;
  for(int i = 0;i<x1.cols();i++)
  {
      Eigen::Vector3d p = x1.col(i);
      double error = n.dot(p) + d;
      error = error < 0 ? (-error):error;
      error_sum += error;
      if(error< maxErrorThreshold)
      {
        num_valid ++;
      }
  }
  return error_sum;
}


double CalculateErrorH(const Eigen::Matrix<double, 2, Eigen::Dynamic>& x1,
                     const Eigen::Matrix<double, 2, Eigen::Dynamic>& x2,
                     const Eigen::Matrix3d& H,const double maxNormThreshold,size_t& num_valid)
{
  double norm_sum = 0.0;
    num_valid = 0;
    for(int i = 0;i<x1.cols();i++)
    {
      Eigen::Vector3d p1(x1.col(i)(0),x1.col(i)(1),1.0);
      Eigen::Vector3d p2(x2.col(i)(0),x2.col(i)(1),1.0);

      Eigen::Vector3d Hp1 = H*p1;
      Hp1 /= Hp1(2);

      double norm_ = (Hp1 - p2).norm();
      if(maxNormThreshold > norm_)
      {
          num_valid ++;
      }
      norm_sum += norm_;

    }
    return norm_sum;
}




}  // namespace colmap
