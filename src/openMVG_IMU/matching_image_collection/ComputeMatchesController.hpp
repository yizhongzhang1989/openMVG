// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_MATCHING_IMAGE_COLLECTION_COMPUTE_MATCHES_CONTROLLER_HPP
#define OPENMVG_IMU_MATCHING_IMAGE_COLLECTION_COMPUTE_MATCHES_CONTROLLER_HPP


#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"

#include <cstdlib>
#include <iostream>

#include <string>


namespace openMVG {
namespace matching_image_collection {
	
	/// Compute corresponding features between a series of views:
	/// - Load view images description (regions: features & descriptors)
	/// - Compute putative local feature matches (descriptors matching)
	/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
	/// - Export computed data
	using namespace sfm;
	class ComputeMatchesController
	{
	public:
		

		static bool Process(const SfM_Data& sfm_data, std::string sMatchesDirectory,
			Matches_Provider* matches_provider, std::string sGeometricModel = "f",
			bool use_predefinedPairs = false, const Pair_Set& sPredefinedPairList = Pair_Set(),
			bool output_matches = false, std::string filepath = "",
			std::string sNearestMatchingMethod = "AUTO", bool bForce = true,
			float fDistRatio = 0.8f, int iMatchingVideoMode = -1,
			bool bGuided_matching = false, int imax_iteration = 2048,
			unsigned int ui_max_cache_size = 0);
	};
}
}

#endif // !1