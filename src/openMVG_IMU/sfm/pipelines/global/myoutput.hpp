// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MY_OUTPUT_HPP
#define MY_OUTPUT_HPP

#include <string>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

bool Output_trajectory(std::string filename,const SfM_Data& sfm_data_);

bool Output_Matchings(std::string filename,const Matches_Provider* matches_provider_);

bool Output_TriangulatedCorrespondings(std::string filename,const SfM_Data& sfm_data_);

bool Output_TriangulatedMatchings(std::string filename,const Matches_Provider* matches_provider_,const SfM_Data& sfm_data_);

bool Output_AngleBetweenRotations(std::string filename,const Matches_Provider* matches_provider_,const SfM_Data& sfm_data_);

void MinMaxMedianMean(std::ofstream& infofile,const std::vector<double>& vec_dif);


} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
