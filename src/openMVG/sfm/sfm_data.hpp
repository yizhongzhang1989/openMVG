// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_HPP
#define OPENMVG_SFM_SFM_DATA_HPP

#include <string>
#include <fstream>
#include <software/SfM/VISfM/VISfM_src/IMU_InteBase.hpp>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"

namespace openMVG {
namespace sfm {

class IMU_Speed
{
public:
    IMU_Speed() = delete;
    IMU_Speed(const Vec3& _speed, double _td):speed_(_speed), td_(_td)
    {
        al_opti = false;
    }
    Vec3 speed_;

    double td_;

    bool al_opti;
};

//using Imus = Hash_Map<IndexT, std::shared_ptr<IMU_InteBase>>;
using Imus = Hash_Map<IndexT, IMU_InteBase>;

using Timestamps = Hash_Map<IndexT, double>;

using SpeedTd = Hash_Map<IndexT, IMU_Speed>;

/// Define a collection of IntrinsicParameter (indexed by View::id_intrinsic)
using Intrinsics = Hash_Map<IndexT, std::shared_ptr<cameras::IntrinsicBase>>;

/// Define a collection of Pose (indexed by View::id_pose)
using Poses = Hash_Map<IndexT, geometry::Pose3>;

using PoseLandmarks = Hash_Map<IndexT, std::set<IndexT>>;

/// Define a collection of View (indexed by View::id_view)
using Views = Hash_Map<IndexT, std::shared_ptr<View>>;

using Tds = Hash_Map<IndexT, double>;

/// Generic SfM data container
/// Store structure and camera properties:
struct SfM_Data
{
  /// Considered views
  Views views;
  /// Considered poses (indexed by view.id_pose)
  Poses poses;

  Poses poses_gt;
  // Considered IMU integration (indexed by view.id_pose)
  Imus imus;

  PoseLandmarks pose_landmark_;

  std::vector<std::pair<double, Pair>> scoring_per_pair_;

  double td_;

  Tds tds;
  // Considered timestamp (indexed by view.id_pose)
  Timestamps timestamps;
  // Considered SpeedBiase (indexed by view.id_pose)
  SpeedTd Speeds;
  std::shared_ptr<IMU_Dataset> imu_dataset;
//  IMU_Dataset imu_dataset;
  Mat3 IG_Ric;
  Vec3 IG_tic;

  std::set<IndexT> add_viewId_cur;


  /// Considered camera intrinsics (indexed by view.id_intrinsic)
  Intrinsics intrinsics;
  /// Structure (3D points with their 2D observations)
  Landmarks structure;
  /// Controls points (stored as Landmarks (id_feat has no meaning here))
  Landmarks control_points;

  /// Root Views path
  std::string s_root_path;

  //--
  // Accessors
  //--
  const Views & GetViews() const {return views;}
  const Poses & GetPoses() const {return poses;}
  const Intrinsics & GetIntrinsics() const {return intrinsics;}
  const Landmarks & GetLandmarks() const {return structure;}
  const Landmarks & GetControl_Points() const {return control_points;}

  /// Check if the View have defined intrinsic and pose
  bool IsPoseAndIntrinsicDefined(const View * view) const
  {
    if (!view) return false;
    return (
      view->id_intrinsic != UndefinedIndexT &&
      view->id_pose != UndefinedIndexT &&
      intrinsics.find(view->id_intrinsic) != intrinsics.end() &&
      poses.find(view->id_pose) != poses.end());
  }

  /// Get the pose associated to a view
  const geometry::Pose3 GetPoseOrDie(const View * view) const
  {
    return poses.at(view->id_pose);
  }
};

struct IMU_Data
{
    Imus imus;
    Timestamps timestamps;
};

struct VISfM_Data : public SfM_Data
{
    Imus imus;
    Timestamps timestamps;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_HPP
