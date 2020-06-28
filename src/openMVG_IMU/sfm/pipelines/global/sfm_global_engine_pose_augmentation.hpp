// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

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
class GlobalSfMReconstructionEngine_PoseAugmentation : public GlobalSfMReconstructionEngine_RelativeMotions_General
{
public:

  GlobalSfMReconstructionEngine_PoseAugmentation(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_PoseAugmentation();

  
  void SetExtraMatchesProvider(Matches_Provider * provider);

  
  
  bool Process();

protected:

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );


  Matches_Provider  * extra_matches_provider_;
  
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
