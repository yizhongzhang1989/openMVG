// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_SFM_GLOBAL_ENGINE_ITERATIVE_POSE_AUGMENTATION_HPP
#define OPENMVG_IMU_SFM_GLOBAL_ENGINE_ITERATIVE_POSE_AUGMENTATION_HPP

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
  void SetMatchesDir(const std::string MatchesDir);

  //iterative augmentation by following
  //	*  loop detection
  //    *  accept the known rotation and refine translations by translation averaging 
  //    *  optimize by ba
  bool Run();

  // accept the known rotation and refine translations by translation averaging 
  // and recompute the structure
  //                          
  // @param[in] extra_matches_provider_    the additional matches
  // the output of function is that update the translation of sfm data.
  // (Taken from OpenMVG with modification)
  bool Process(std::shared_ptr<Matches_Provider> extra_matches_provider_);


  // find the proper camera pair whose angular is less than 25 degree.
  //                          
  // @param[out] extra_pairs       container of output paris
  // @param[in]  tried_pairs       the pairs tried once
  // (Owned by BC)
  bool LoopDetection(Pair_Set& extra_pairs, const Pair_Set& tried_pairs);
  bool Optimize();
protected:

  // find the proper camera pair whose angular is less than 25 degree.
  //                          
  // @param[in]  global_rotations          global rotation of every image
  // @param[out] tripletWise_matches       the inlier matches that comfirm to the translations.
  // @param[in] extra_matches_provider_    the additional matches
  // the output of function is that update the translation of sfm data.
  // (Taken from OpenMVG with modification)
  bool Compute_Global_PrioTranslations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches,
    Matches_Provider* extra_matches_provider_
  );

  std::string sMatchesDir_;

  
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
