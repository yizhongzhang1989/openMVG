// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
/// Spurious correspondences are discarded by using the
///  a threshold over the distance ratio of the 2 nearest neighbours.
/// Using a Cascade Hashing matching
/// Cascade hashing tables are computed once and used for all the regions.
///
class Hierarchical_Matcher_Regions : public Matcher
{
  public:

	explicit Hierarchical_Matcher_Regions
	(
		float distRatio, std::string bin_dir, double maxDistanceThreshold, const sfm::SfM_Data& sfm_data,
		const std::string sMatches_dir, bool bfeature_validation = true, bool bopticalfiltering = true,
		bool bdynamicdistance = true, bool bopticalmatching = true, bool bdebug = true
	);

  /// Find corresponding points between some pair of view Ids
  void Match
  (const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair_Set & pairs,
    matching::PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
    C_Progress * progress = nullptr
  ) const override;

  void Hierarchical_Match
  (const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
	  const Pair_Set & pairs,
	  matching::PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
	  C_Progress * progress = nullptr
  );

  



 public:
  // Distance ratio used to discard spurious correspondence
  float f_dist_ratio_;
  std::string sMatches_dir_;  // for debug
  matching::OpticalFlow_Container opticalflow_container;
  std::map<IndexT, std::set<IndexT>> notable_features;
  bool bfeature_validation_;
  bool bopticalfiltering_;
  bool bdynamicdistance_;
  bool bopticalmatching_;
  bool bdebug_;
  
};

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_CASCADE_HASHING_MATCHER_REGIONS_HPP
