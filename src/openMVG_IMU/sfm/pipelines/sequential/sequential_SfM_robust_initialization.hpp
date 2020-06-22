// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_ROBUST_INITIALIZATION_HPP
#define OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_ROBUST_INITIALIZATION_HPP

#include <set>
#include <string>
#include <vector>

#include "openMVG_IMU/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/multiview/triangulation_method.hpp"
#include "openMVG/tracks/tracks.hpp"

//namespace htmlDocument { class htmlDocumentStream; }
//namespace { template <typename T> class Histogram; }

namespace openMVG {
namespace sfm {

struct Features_Provider;
struct Matches_Provider;

/// Sequential SfM Pipeline Reconstruction Engine.
class SequentialSfMReconstructionEngine_Robust_Initialization : public SequentialSfMReconstructionEngine_General
{
public:

  SequentialSfMReconstructionEngine_Robust_Initialization(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "",
    const size_t initial_max_iteration_count = 256);

  ~SequentialSfMReconstructionEngine_Robust_Initialization();

  bool Process_Robust_Initialization();
  bool Process_KnownTracks();

  bool RobustAutomaticInitialPairChoice(Pair & initial_pair) const;

  /// Return MSE (Mean Square Error) and a histogram of residual values.
  double ComputeResidualsHistogram(Histogram<double> * histo);

  /// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
  bool RobustMakeInitialPair3D(const Pair & initialPair);

  void setIMUData(const SfM_Data & imu_data)
  {
    imu_data_ = imu_data;
    b_robust_initialization_of_imu_ = true;
  }
  

protected:


public:    //bc
    size_t initial_max_iteration_count_;
    SfM_Data imu_data_;
    bool b_robust_initialization_of_imu_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_HPP
