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
namespace utils{
using namespace sfm;

bool Output_trajectory(std::string filename,const SfM_Data& sfm_data_);

bool Output_Matches_provider(std::string filename,const Matches_Provider* matches_provider_);

bool Output_TriangulatedCorrespondings(std::string filename,const SfM_Data& sfm_data_);

bool Output_TriangulatedMatchings(std::string filename,const Matches_Provider* matches_provider_,const SfM_Data& sfm_data_);

bool Output_AngleBetweenRotations(std::string filename,const Matches_Provider* matches_provider_,const SfM_Data& sfm_data_);


template<typename T>
void MinMaxMedianMean(std::ostream& infofile,const std::vector<T>& vec_dif)
{
    if(vec_dif.size()==0) 
    {
        infofile<<"Warning:the vector is empty\n";
        return;
    }
    T sum = vec_dif[0];
    for(size_t i = 1 ;i <vec_dif.size(); i++ )
    {
        sum += vec_dif[i];
        if(vec_dif[i]<vec_dif[i-1])
        {
            infofile<<"Error:the difference vector is not sorted\n";
            return;
        }
    }
    infofile<<"Min:"<<vec_dif[0]<<"\n";
    infofile<<"Max:"<<*(vec_dif.rbegin())<<"\n";
    infofile<<"Median:"<<vec_dif[vec_dif.size()/2]<<"\n";
    infofile<<"Mean:"<<sum/((double)vec_dif.size())<<"\n";
	infofile<< "Sum:" << sum << "\n";

}

template<typename T>
void AllMinMaxMedianMean(std::ostream& infofile, const std::vector<T>& vec_dif)
{
	if (vec_dif.size() == 0)
	{
		infofile << "Warning:the vector is empty\n";
		return;
	}
	T sum = vec_dif[0];
	infofile << vec_dif[0] << "\n";
	for (size_t i = 1; i <vec_dif.size(); i++)
	{
		sum += vec_dif[i];
		infofile << vec_dif[i] << "\n";
		if (vec_dif[i]<vec_dif[i - 1])
		{
			infofile << "Error:the difference vector is not sorted\n";
			return;
		}
	}
	infofile << "Min:" << vec_dif[0] << "\n";
	infofile << "Max:" << *(vec_dif.rbegin()) << "\n";
	infofile << "Median:" << vec_dif[vec_dif.size() / 2] << "\n";
	infofile << "Mean:" << sum / ((double)vec_dif.size()) << "\n";
	infofile << "Sum:" << sum << "\n";

}

inline bool Output_Matches(std::ostream& os, const matching::PairWiseMatches& map_PutativeMatches, bool dif_frames_analysis = false)
{
	os << "Output_Matches\n";
	std::map<IndexT, std::vector<IndexT>> statistic_pairs;
	for (const auto& pairwise_match : map_PutativeMatches)
	{
		os << pairwise_match.first.first << "," << pairwise_match.first.second << ":" << pairwise_match.second.size() << "\n";
		IndexT dif_frames = pairwise_match.first.second > pairwise_match.first.first ?
			(pairwise_match.first.second - pairwise_match.first.first) :
			(pairwise_match.first.first - pairwise_match.first.second);

		statistic_pairs[dif_frames].emplace_back(pairwise_match.second.size());
	}
	if (dif_frames_analysis)
	{
		for (auto& statistic_pair_item : statistic_pairs)
		{
			os << "statistic of matches with " << statistic_pair_item.first << " apart frames\n";
			std::sort(statistic_pair_item.second.begin(), statistic_pair_item.second.end());
			MinMaxMedianMean<IndexT>(os, statistic_pair_item.second);
		}
	}

	return true;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
