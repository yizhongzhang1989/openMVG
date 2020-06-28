// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
#define OPENMVG_IMU_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP

#include <string>
#include <vector>

#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"  //BC
#include "openMVG/tracks/tracks.hpp"

namespace openMVG { namespace graph { struct Triplet; } }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{


struct SfM_Data;
struct Matches_Provider;
struct Features_Provider;
/// Note:the class is created as same as the homonymous class in openMVG,
///      but the only difference is that all member variables and functions
///      in the class are declared as public for inheriting.
class GlobalSfM_Translation_AveragingSolver_General
{
  

public:
	std::vector<RelativeInfo_Vec> vec_relative_motion_;
	////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  const std::vector<RelativeInfo_Vec> & Getrelative_motion() const {return vec_relative_motion_;} 
  //END(Author: BC)===================================================

  bool Run(
	  ETranslationAveragingMethod eTranslationAveragingMethod,
    openMVG::sfm::SfM_Data & sfm_data,
    const openMVG::sfm::Features_Provider * features_provider,
    const openMVG::sfm::Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches & tripletWise_matches
  );

public:   //BC
  bool Translation_averaging(
	  ETranslationAveragingMethod eTranslationAveragingMethod,
    sfm::SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR);

  void Compute_translations(
    const sfm::SfM_Data & sfm_data,
    const sfm::Features_Provider * features_provider,
    const sfm::Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches &tripletWise_matches);

  //-- Compute the relative translations on the rotations graph.
  // Compute relative translations by using triplets of poses.
  // Use an edge coverage algorithm to reduce the graph covering complexity
  // Complexity: sub-linear in term of edges count.
  void ComputePutativeTranslation_EdgesCoverage(
    const sfm::SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const sfm::Features_Provider * features_provider,
    const sfm::Matches_Provider * matches_provider,
    std::vector<RelativeInfo_Vec> & vec_triplet_relative_motion,
    matching::PairWiseMatches & newpairMatches);

  // Robust estimation and refinement of triplet of translations
  bool Estimate_T_triplet(
    const sfm::SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const sfm::Features_Provider * features_provider,
    const sfm::Matches_Provider * matches_provider,
    const graph::Triplet & poses_id,
    std::vector<Vec3> & vec_tis,
    double & dPrecision, // UpperBound of the precision found by the AContrario estimator
    std::vector<uint32_t> & vec_inliers,
    openMVG::tracks::STLMAPTracks & rig_tracks,
    const std::string & sOutDirectory) const;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
