// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMU_SFM_GLOBAL_ENGINE_POSE_AUGMENTATION_HPP
#define OPENMVG_IMU_SFM_GLOBAL_ENGINE_POSE_AUGMENTATION_HPP

#include <memory>
#include <string>

#include "openMVG_IMU/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG_IMU/sfm/pipelines/global/GlobalSfM_priotranslation_averaging.hpp"
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
class GlobalSfMReconstructionEngine_PoseAugmentation : public ReconstructionEngine
{
public:

  GlobalSfMReconstructionEngine_PoseAugmentation(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_PoseAugmentation() override;

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);
  void SetExtraMatchesProvider(Matches_Provider * provider);

  void SetRotationAveragingMethod(ERotationAveragingMethod eRotationAveragingMethod);
  void SetTranslationAveragingMethod(ETranslationAveragingMethod eTranslation_averaging_method_);
  
  bool Process() override;

protected:

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );

  /// Compute the initial structure of the scene
  bool Compute_Initial_Structure
  (
    matching::PairWiseMatches & tripletWise_matches
  );




  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  ERotationAveragingMethod eRotation_averaging_method_;
  ETranslationAveragingMethod eTranslation_averaging_method_;

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;
  Matches_Provider  * extra_matches_provider_;
  
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
