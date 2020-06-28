// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_SFM_GLOBAL_ENGINE_MOTION_AVERAGING_HPP
#define OPENMVG_IMU_SFM_GLOBAL_ENGINE_MOTION_AVERAGING_HPP

#include <memory>
#include <string>


#include "openMVG_IMU/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG_IMU/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_MotionAveraging : public GlobalSfMReconstructionEngine_RelativeMotions_General
{
public:

  GlobalSfMReconstructionEngine_MotionAveraging(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_MotionAveraging() ;

  bool Process() ;

protected:
  /// Compute from relative rotations the global rotations of the camera poses

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations_PrintMotion
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );

  /// Compute the initial structure of the scene
  

  // Adjust the scene (& remove outliers)


private:
  /// Compute relative rotations
  

  //----
  //-- Data
  //----

};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
