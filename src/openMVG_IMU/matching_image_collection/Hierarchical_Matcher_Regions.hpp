// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_MATCHING_HIERARCHICAL_MATCHER_REGIONS_HPP
#define OPENMVG_IMU_MATCHING_HIERARCHICAL_MATCHER_REGIONS_HPP

#include <memory>

#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG_IMU/matching/optical_flow.hpp"

namespace openMVG { namespace matching { class PairWiseMatchesContainer; } }
namespace openMVG { namespace sfm { struct Regions_Provider; } }

namespace openMVG {
namespace matching_image_collection {

/// Implementation of an Hierarchical Matcher
/// sift matching (+feature validation) (+static/dynamic optical filtering)
///                   (+static/dynamic optical matching).
///
class Hierarchical_Matcher_Regions : public Matcher
{
  public:

	explicit Hierarchical_Matcher_Regions
	(
		float distRatio, std::string bin_dir, double maxDistanceThreshold, const sfm::SfM_Data& sfm_data,
		bool bfeature_validation = false, bool bopticalfiltering = true,
		bool bdynamicdistance = false, bool bopticalmatching = false, 
		bool bnotablefeaturesdetection = false, bool bnotablevalidation = false
	);

  /// Find corresponding points between some pair of view Ids
  void Match
  (const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair_Set & pairs,
    matching::PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
    C_Progress * progress = nullptr
  ) const override;

  void Match_Hierarchical
  (
	  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
	  const Pair_Set & pairs,
	  matching::PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
	  C_Progress * my_progress_bar
  );

  



 public:
  // Distance ratio used to discard spurious correspondence
  float f_dist_ratio_;
  
   
  // Process optical informations
  matching::OpticalFlow_Container opticalflow_container;
  // notable features
  std::map<IndexT, std::set<IndexT>> notable_features;
  bool bfeature_validation_;
  bool bopticalfiltering_;
  bool bdynamicdistance_;
  bool bopticalmatching_;
  bool bnotablefeaturesdetection_;
  bool bnotablevalidation_;
  
};

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_CASCADE_HASHING_MATCHER_REGIONS_HPP
