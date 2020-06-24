// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_MATCHING_OPTICAL_FLOW_MATCHER_REGIONS_HPP
#define OPENMVG_IMU_MATCHING_OPTICAL_FLOW_MATCHER_REGIONS_HPP

#include <memory>

#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG_IMU/matching/optical_flow.hpp"

namespace openMVG { namespace matching { class PairWiseMatchesContainer; } }
namespace openMVG { namespace sfm { struct Regions_Provider; } }

namespace openMVG {
namespace matching_image_collection {

/// Implementation of an Optical Matcher
/// Compute putative matches between a collection of pictures
/// Spurious correspondences are discarded by using the
///  a threshold over the distance ratio of the 2 nearest neighbours.
/// Using a Cascade Hashing matching
/// Cascade hashing tables are computed once and used for all the regions.
///(Owned by BC)
class Optical_Flow_Matcher_Regions : public Matcher
{
  public:

	struct Optical_track
	{
		unsigned int view_id;
		unsigned int feat_id;
		float x;
		float y;
		float rx;
		float ry;
	};
	explicit Optical_Flow_Matcher_Regions
	(
		float distRatio, std::string bin_dir, double maxDistanceThreshold, const sfm::SfM_Data& sfm_data,
		std::string sMatches_dir
	);

  /// Find corresponding points between some pair of view Ids
  void Match
  (
	  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
	  const Pair_Set & pairs,
	  matching::PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
      C_Progress * progress = nullptr
  ) const override;

  private:
  // Distance ratio used to discard spurious correspondence
  float f_dist_ratio_;
  // the directory stores the binary file of features tracked in optical flow.
  std::string bin_dir_;
  // the distance threshold for optical filtering / matching
  double MaxDistanceThreshold;
  std::string sMatches_dir_;
  // Process optical informations
  matching::OpticalFlow_Container opticalflow_container;
};

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_CASCADE_HASHING_MATCHER_REGIONS_HPP
