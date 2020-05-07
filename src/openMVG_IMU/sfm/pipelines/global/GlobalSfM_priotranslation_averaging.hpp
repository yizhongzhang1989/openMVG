// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMU_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_PRIOTRANSLATION_AVERAGING_HPP  //BC
#define OPENMVG_IMU_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_PRIOTRANSLATION_AVERAGING_HPP

#include <string>
#include <vector>


#include "openMVG_IMU/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"  //BC


namespace openMVG { namespace graph { struct Triplet; } }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{


//BC START//
class GlobalSfM_PrioTranslation_AveragingSolver: public GlobalSfM_Translation_AveragingSolver_General 
{

public:
  bool MyRun
  (
    ETranslationAveragingMethod eTranslationAveragingMethod,
    openMVG::sfm::SfM_Data & sfm_data,
    const openMVG::sfm::Features_Provider * features_provider,
    const openMVG::sfm::Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches & tripletWise_matches,
    const openMVG::sfm::Matches_Provider * extra_matches_provider
  );

public:   //BC
  

  bool PrioTranslation_averaging(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  sfm::SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR);

  //-- Compute the relative translations on the rotations graph.
  // Compute relative translations by using triplets of poses.
  // Use an edge coverage algorithm to reduce the graph covering complexity
  // Complexity: sub-linear in term of edges count.
  
};
//BC END//
} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
