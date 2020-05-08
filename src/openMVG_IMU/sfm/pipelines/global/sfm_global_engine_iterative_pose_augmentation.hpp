// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMU_SFM_GLOBAL_ENGINE_POSE_AUGMENTATION_HPP
#define OPENMVG_IMU_SFM_GLOBAL_ENGINE_POSE_AUGMENTATION_HPP

#include <memory>
#include <string>

#include "openMVG_IMU/sfm/pipelines/global/GlobalSfM_priotranslation_averaging.hpp"
#include "openMVG_IMU/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"  //BC
#include "openMVG/sfm/pipelines/sfm_engine.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_IterativePoseAugmentation : public GlobalSfMReconstructionEngine_RelativeMotions_General
{
public:

  GlobalSfMReconstructionEngine_IterativePoseAugmentation(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_IterativePoseAugmentation();

  
  bool Run();
  bool Process();
  bool LoopDetection();
protected:

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );

  

  
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
