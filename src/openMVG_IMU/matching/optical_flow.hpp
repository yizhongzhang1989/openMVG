// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#ifndef OPENMVG_IMU_MATCHING_OPTICAL_FILTERING_HPP
#define OPENMVG_IMU_MATCHING_OPTICAL_FILTERING_HPP




#include <cmath>
#include <random>
#include <utility>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "openMVG/types.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/cascade_hasher.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG {
namespace matching {

//Optical Information Controller Class
//    *	optical filtering
//    * optical matching
//    * dynamic selection of threshold

class OpticalFlow_Container {

public:
//The data format in binary file of optical flow
struct Optical_track
{
	unsigned int view_id;
	unsigned int feat_id;
	float x;
	float y;
	float rx;
	float ry;
};

double double_comp(const double x, const double y) const
{
	const double eps = 1e-5;
	if (x - y >= -eps&&x - y <= eps)
	{
		return 0;
	}
	else if (x - y > eps)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

double cal_dis (double x1, double y1, double x2, double y2) const
{
		return std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
};

explicit OpticalFlow_Container(const std::string bin_dir,
								const double MaxDistanceThreshold,
								const sfm::SfM_Data& sfm_data)
	 :bin_dir_(bin_dir),MaxDistanceThreshold_(MaxDistanceThreshold)
  {
	  readable_ = false;
	//read binary data of optical flow
	std::cout << "**read binary features**\n";
	for (const auto& view_item : sfm_data.GetViews())
	{
		const IndexT view_id = view_item.first;
		std::stringstream ss;
		Optical_track ot;
		//read the binary file by view_id
		//ss << std::setw(5) << std::setfill('0') << std::stoi(stlplus::basename_part(view_item.second->s_Img_path)) << ".bin";
		ss << std::setw(5) << std::setfill('0') << view_id << ".bin";
		std::string filename = ss.str();
		std::ifstream file(stlplus::create_filespec(bin_dir, filename), std::ios::binary);

		if (!file.is_open())
		{
			std::cerr << stlplus::create_filespec(bin_dir, std::string(filename)) << " can not be read\n";
			return ;
		}
		if (opticaltrack_table_.count(view_id))
		{
			std::cerr << "Error:duplicate view id\n";
			return ;
		}
		opticaltrack_table_.emplace(view_id, std::map<Pair, Optical_track>());
		std::map<Pair, Optical_track>& map_optrack = opticaltrack_table_.at(view_id);
		// read the position information from file
		while (file.read(reinterpret_cast<char*>(&ot), sizeof(Optical_track)))
		{
			map_optrack[Pair(ot.view_id, ot.feat_id)] = ot;
			if (file.peek() == '\n')
			{
				file.ignore();
			}
		}
		
		std::cout << "#features of view  " << view_id << "("<< filename <<"): " << opticaltrack_table_.at(view_id).size() << "\n";

	}
	  readable_ = true;
  }

// compute matches in local and for every feature in first view ,find least NN features paired in second view.
//      view1											view2
//    -------------------						-------------------
//	  | feat_id1(view1)-+-----------------------+-feat_id2(view2) |
//    |        |        |    feature match		|       |---------+---> distance tracked by optical flow
//    |        ---------+-----------------------+>feat_id1(view1) |
//    -------------------	optical flow		-------------------
//                          
//
// @param[in]  view_id_1                      the first view id
// @param[in]  hashed_descriptions1           the hased descriptors of first view
// @param[in]  descriptions1                  the descriptors of first view(N*M,N:number of features;M:length of descriptors)
// @param[in]  view_id_2                      the second view id
// @param[in]  hashed_descriptions2           the hased descriptors of second view
// @param[in]  descriptions2                  the descriptors of second view(N*M,N:number of features;M:length of descriptors)
// @param[in]  pts_1                          the feature positions of first view(used for validating correctness)
// @param[in]  pts_2                          the feature positions of first view(used for validating correctness)
// @param[in]  latter_view_id                 the back view id among two views.
// @param[in]  nb_hash_code                   the bits of hash code
// @param[out] pvec_indices                   the matched feature pair
// @param[out] pvec_distances                 distance of every matches feature pair
// @param[in]  paired_feati                   the mask of features in first view
// @param[in]  paired_featj                   the mask of features in second view
// @param[in]  NN                             the required number of features paired for one feature
// (Taken from openMVG with modification)

template <typename MatrixT, typename DistanceType>
bool Optical_Matching
(
	const IndexT view_id_1,
	const HashedDescriptions& hashed_descriptions1,
	const MatrixT & descriptions1,
	const IndexT view_id_2,
	const HashedDescriptions& hashed_descriptions2,
	const MatrixT & descriptions2,
	const std::vector<features::PointFeature>& pts_1,
	const std::vector<features::PointFeature>& pts_2,
	
	const IndexT latter_view_id,
	const uint8_t nb_hash_code,
	IndMatches * pvec_indices,
	std::vector<DistanceType> * pvec_distances,
	
	const int NN = 2,
	const std::set<IndexT>* paired_feati = nullptr,
	const std::set<IndexT>* paired_featj = nullptr
)const
{

	const std::map<Pair, Optical_track>& latter_opticaltable = opticaltrack_table_.at(latter_view_id);
	using MetricT = L2<typename MatrixT::Scalar>;
	MetricT metric;
	
	const int kNumTopCandidates = 10;
	

	//optical flow match
	//Is the feature id consistent with index of vector<features::PointFeature> and column of descriptor matrix
	//

	// Preallocated hamming distances. Each column indicates the hamming distance
	// and the rows collect the descriptor ids with that
	// distance. num_descriptors_with_hamming_distance keeps track of how many
	// descriptors have that distance.
	Eigen::MatrixXi candidate_hamming_distances(
			hashed_descriptions2.hashed_desc.size(), nb_hash_code + 1);
	
	Eigen::VectorXi num_descriptors_with_hamming_distance(nb_hash_code + 1);
	// Preallocate the container for keeping euclidean distances.
	std::vector<std::pair<DistanceType, int>> candidate_euclidean_distances;
	candidate_euclidean_distances.reserve(kNumTopCandidates);

	using HammingMetricType = matching::Hamming<stl::dynamic_bitset::BlockType>;
	static const HammingMetricType metricH = {};
	
	for (size_t i = 0; i < hashed_descriptions1.hashed_desc.size(); i++)
	{

		if(paired_feati != nullptr && paired_feati->count(i)) continue;
		num_descriptors_with_hamming_distance.setZero();
		candidate_euclidean_distances.clear();
		////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
		size_t feats_selected = 0;
		//check whether the first feature in first view is tracked.
		if (!latter_opticaltable.count(Pair(view_id_1, i)))
		{
			if (latter_view_id == view_id_1)
			{
				//the first feature  is not tracked in first view
				std::cout << "Error:feature " << i << "can not be found in the binary file of view " << view_id_1 << "\n";
				return false;
			}
			else
			{
				//the first feature  is not tracked in second view
				continue;
			}
		}
		const Optical_track& ot_1 = latter_opticaltable.at(Pair(view_id_1, i));
		const auto& hashed_desc = hashed_descriptions1.hashed_desc[i];
		
		for (size_t j = 0; j < hashed_descriptions2.hashed_desc.size(); j++)
		{
			if(paired_featj !=nullptr && paired_featj->count(j)) continue;
			//check whether the second feature in second view is tracked.
			if (!latter_opticaltable.count(Pair(view_id_2, j)))
			{
				if (latter_view_id == view_id_2)
				{
					//the second feature  is not tracked in second view
					std::cout << "Error:feature " << j << "can not be found in the binary file of view " << view_id_2 << "\n";
					return false;
				}
				else
				{
					//the second feature  is not tracked in first view
					continue;
				}
			}
			//check the correctness of features
			const Optical_track& ot_2 = latter_opticaltable.at(Pair(view_id_2, j));
			if (double_comp(ot_1.rx, pts_1[i].x()) != 0 ||
				double_comp(ot_1.ry, pts_1[i].y()) != 0 ||
				double_comp(ot_2.rx, pts_2[j].x()) != 0 ||
				double_comp(ot_2.ry, pts_2[j].y()) != 0)
			{
				std::cerr << "Error: optical features is inconsistent with original features\n";
				return false;
			}
			// compute the tracked distance of two features
			double dis_featpair = cal_dis(ot_1.x, ot_1.y, ot_2.x, ot_2.y);
			// compute hamming distance of feature pair whose optical distance is less than `MaxDistanceThreshold_`
			// and Put the descriptors into buckets corresponding to their hamming distance.
			if (dis_featpair <= MaxDistanceThreshold_)
			{
				const HammingMetricType::ResultType hamming_distance = metricH(
					hashed_desc.hash_code.data(),
					hashed_descriptions2.hashed_desc[j].hash_code.data(),
					hashed_desc.hash_code.num_blocks());
				candidate_hamming_distances(
					num_descriptors_with_hamming_distance(hamming_distance)++,
					hamming_distance) = j;
				feats_selected++;
			}

		}
		//END(Author: BC)===================================================

		// Compute the euclidean distance of the k descriptors with the best hamming
		// distance 
		candidate_euclidean_distances.reserve(kNumTopCandidates);
		for (int j = 0; j < candidate_hamming_distances.cols() &&
			(candidate_euclidean_distances.size() < kNumTopCandidates); ++j)
		{
			for (int k = 0; k < num_descriptors_with_hamming_distance(j) &&
				(candidate_euclidean_distances.size() < kNumTopCandidates); ++k)
			{
				const int candidate_id = candidate_hamming_distances(k, j);
				const DistanceType distance = metric(
					descriptions2.row(candidate_id).data(),
					descriptions1.row(i).data(),
					descriptions1.cols());

				candidate_euclidean_distances.emplace_back(distance, candidate_id);
			}
		}
		// Assert that each query is having at least NN retrieved neighbors
		if (candidate_euclidean_distances.size() >= NN)
		{
			
			// Find the top NN candidates based on euclidean distance.
			std::partial_sort(candidate_euclidean_distances.begin(),
				candidate_euclidean_distances.begin() + NN,
				candidate_euclidean_distances.end());
			// save resulting neighbors
			for (int l = 0; l < NN; ++l)
			{
				pvec_distances->emplace_back(candidate_euclidean_distances[l].first);
				pvec_indices->emplace_back(IndMatch(i, candidate_euclidean_distances[l].second));
			}
		}
	}
	
	return true;

}

// Filtering the matches not satisfying the distance threshold in optical flow between two views
//      view1											view2
//    -------------------						-------------------
//	  | feat_id1(view1)-+-----------------------+-feat_id2(view2) |
//    |        |        |    feature match		|       |---------+---> distance tracked by optical flow
//    |        ---------+-----------------------+>feat_id1(view1) |
//    -------------------	optical flow		-------------------
//                          
//
// @param[in]  first_view_id              the first view id
// @param[in]  second_view_id             the second view id
// @param[out] vec_putative_matches       putative matches to be filtered
// @param[in]  pts_1                      the feature positions of first view(used for validating correctness)
// @param[in]  pts_2                      the feature positions of first view(used for validating correctness)
// @param[in]  latter_view_id             the back view id among two views.
// @param[out] dis_featpairs              the optical distances of all feature pairs
// (Owned by BC)

bool Optical_Filtering
(
	const IndexT first_view_id,
	
	const IndexT second_view_id,
	matching::IndMatches& vec_putative_matches,
	const std::vector<features::PointFeature>& pts_1,
	const std::vector<features::PointFeature>& pts_2,
	const IndexT latter_view_id,

	std::vector<double>* dis_featpairs = nullptr   
)const
{
	// get the tracked feature information in back view.
	const std::map<Pair, Optical_track>& latter_opticaltable = opticaltrack_table_.at(latter_view_id);
	IndMatches::iterator im_iter;
	for (im_iter = vec_putative_matches.begin(); im_iter != vec_putative_matches.end();)
	{
		const IndexT first_feat_id = im_iter->i_;
		const IndexT second_feat_id = im_iter->j_;
		Vec2f first_feat = pts_1[first_feat_id].coords();
		Vec2f second_feat = pts_2[second_feat_id].coords();
		//check whether the first feature in first view is tracked.
		if (latter_opticaltable.count(Pair(first_view_id, first_feat_id)) == 0)
		{
			if (latter_view_id == first_view_id)
			{
				//the first feature  is not tracked in first view
				std::cerr << "Error: feature" << first_view_id << " is losed in view " << first_view_id << "\n";
				return false;
			}
			else
			{
				//the first feature  is not tracked in second view
				im_iter = vec_putative_matches.erase(im_iter);
				
				continue;
			}

		}
		//check whether the second feature in second view is tracked.
		if (latter_opticaltable.count(Pair(second_view_id, second_feat_id)) == 0)
		{
			if (latter_view_id == second_view_id)
			{
				//the second feature  is not tracked in second view
				std::cerr << "Error: feature" << second_feat_id << " is losed in view " << second_view_id << "\n";
				return false;
			}
			else
			{
				//the second feature  is not tracked in first view
				im_iter = vec_putative_matches.erase(im_iter);
				
				continue;
			}
		}
		//check the correctness of features
		const Optical_track& ot_1 = latter_opticaltable.at(Pair(first_view_id, first_feat_id));
		const Optical_track& ot_2 = latter_opticaltable.at(Pair(second_view_id, second_feat_id));

		if (double_comp(ot_1.rx, first_feat[0]) != 0 ||
			double_comp(ot_1.ry, first_feat[1]) != 0 ||
			double_comp(ot_2.rx, second_feat[0]) != 0 ||
			double_comp(ot_2.ry, second_feat[1]) != 0)
		{
			std::cerr << "Error: optical features is inconsistent with raw features\n";
			return false;
		}
		// compute the tracked distance of two features
		double dis_featpair = cal_dis(ot_1.x, ot_1.y, ot_2.x, ot_2.y);
		//record the distance 
		if(dis_featpairs!=nullptr) dis_featpairs->push_back(dis_featpair);

		if (dis_featpair > MaxDistanceThreshold_)   //delete the feature pairs whose distance is too large.
		{
			im_iter = vec_putative_matches.erase(im_iter);
		}
		else
		{
			im_iter++;
		}


	}
	return true;
}

// Filtering the matches not satisfying the distance threshold in optical flow 
//      view1											view2
//    -------------------						-------------------
//	  | feat_id1(view1)-+-----------------------+-feat_id2(view2) |
//    |        |        |    feature match		|       |---------+---> distance tracked by optical flow
//    |        ---------+-----------------------+>feat_id1(view1) |
//    -------------------	optical flow		-------------------
//                          
//
// @param[out] map_PutativesMatches      putative matches to be filtered
// @param[in]  regions_provider          the information of features.
// @param[out] my_progress_bar           the UI of progress bar
// @param[out] dis_featpairs             the optical distances of all feature pairs

// (Owned by BC)

bool Optical_Filtering(matching::PairWiseMatches& map_PutativesMatches,
	const sfm::Regions_Provider* regions_provider, 
	C_Progress * my_progress_bar = nullptr,
	std::vector<double>* dis_featpairs = nullptr   // for debugging
	)const
{
	std::cout << "**filter putative matches**\n";
	if (!my_progress_bar)
		my_progress_bar = &C_Progress::dummy();
	my_progress_bar->restart(map_PutativesMatches.size(), "\n- Optical Filtering -\n");
	

	
	PairWiseMatches::iterator pwm_iter;
	for (pwm_iter = map_PutativesMatches.begin(); pwm_iter != map_PutativesMatches.end(); pwm_iter++)
	{
		if (my_progress_bar->hasBeenCanceled())
			break;
		const IndexT first_view_id = pwm_iter->first.first;
		const IndexT second_view_id = pwm_iter->first.second;
		const IndexT latter_view_id = std::max(first_view_id, second_view_id);
		std::shared_ptr<features::Regions> feats_1 = regions_provider->get(first_view_id);
		std::shared_ptr<features::Regions> feats_2 = regions_provider->get(second_view_id);

		//std::cout << "Pair " << first_view_id << "-" << second_view_id << "\n";
		
		IndMatches::iterator im_iter;
		// get the tracked feature information in back view.
		const std::map<Pair, Optical_track>& latter_opticaltable = opticaltrack_table_.at(latter_view_id);
		// filter the feature matches whose distance is larger than `MaxDistanceThreshold`
		//std::cout << "	#raw feature pairs" << pwm_iter->second.size() << "\n";
		for (im_iter = pwm_iter->second.begin(); im_iter != pwm_iter->second.end();)
		{
			if (my_progress_bar->hasBeenCanceled())
				continue;
			const IndexT first_feat_id = im_iter->i_;
			const IndexT second_feat_id = im_iter->j_;
			Vec2 first_feat = feats_1->GetRegionPosition(first_feat_id);
			Vec2 second_feat = feats_2->GetRegionPosition(second_feat_id);
			//check whether the first feature in first view is tracked.
			if (latter_opticaltable.count(Pair(first_view_id, first_feat_id)) == 0)
			{
				if (latter_view_id == first_view_id)
				{
					//the first feature  is not tracked in first view
					std::cerr << "Error: feature" << first_view_id << " is losed in view " << first_view_id << "\n";
					return false;
				}
				else
				{
					//the first feature  is not tracked in second view
					im_iter = pwm_iter->second.erase(im_iter);
					
					continue;
				}

			}
			//check whether the second feature in second view is tracked.
			if (latter_opticaltable.count(Pair(second_view_id, second_feat_id)) == 0)
			{
				if (latter_view_id == second_view_id)
				{
					//the second feature  is not tracked in second view
					std::cerr << "Error: feature" << second_feat_id << " is losed in view " << second_view_id << "\n";
					return false;
				}
				else
				{
					//the second feature  is not tracked in first view
					im_iter = pwm_iter->second.erase(im_iter);
					
					continue;
				}
			}
			//check the correctness of features
			const Optical_track& ot_1 = latter_opticaltable.at(Pair(first_view_id, first_feat_id));
			const Optical_track& ot_2 = latter_opticaltable.at(Pair(second_view_id, second_feat_id));

			if (double_comp(ot_1.rx, first_feat[0]) != 0 ||
				double_comp(ot_1.ry, first_feat[1]) != 0 ||
				double_comp(ot_2.rx, second_feat[0]) != 0 ||
				double_comp(ot_2.ry, second_feat[1]) != 0)
			{
				std::cerr << "Error: optical features is inconsistent with raw features\n";
				return false;
			}
			// compute the tracked distance of two features
			double dis_featpair = cal_dis(ot_1.x, ot_1.y, ot_2.x, ot_2.y);
			if(dis_featpairs!=nullptr) dis_featpairs->push_back(dis_featpair);
			if (dis_featpair > MaxDistanceThreshold_)   //delete the feature pairs whose distance is too large.
			{
				im_iter = pwm_iter->second.erase(im_iter);
			}
			else
			{
				im_iter++;
			}


		}

		++(*my_progress_bar);

	}
	
	return true;
}

// Compute the distance of all feature pairs to set a dynamic distance threshold
//
// @param[in] map_PutativesMatches      putative matches to be filtered
// @param[in]  regions_provider          the information of features.


// (Owned by BC)

bool SetDynamicDistanceThreshold(const matching::PairWiseMatches& map_PutativesMatches,
	const sfm::Regions_Provider* regions_provider
)
{
	// record the distance informations.
	std::vector<double> dis_featpairs;

	std::cout<<"calculating distance threshold\n";
	PairWiseMatches::const_iterator pwm_iter;
	for (pwm_iter = map_PutativesMatches.begin(); pwm_iter != map_PutativesMatches.end(); pwm_iter++)
	{
		const IndexT first_view_id = pwm_iter->first.first;
		const IndexT second_view_id = pwm_iter->first.second;
		const IndexT latter_view_id = std::max(first_view_id, second_view_id);
		std::shared_ptr<features::Regions> feats_1 = regions_provider->get(first_view_id);
		std::shared_ptr<features::Regions> feats_2 = regions_provider->get(second_view_id);

		//std::cout << "Pair " << first_view_id << "-" << second_view_id << "\n";

		IndMatches::const_iterator im_iter;
		// get the tracked feature information in back view.
		const std::map<Pair, Optical_track>& latter_opticaltable = opticaltrack_table_.at(latter_view_id);
	
		for (im_iter = pwm_iter->second.begin(); im_iter != pwm_iter->second.end(); im_iter++)
		{
			const IndexT first_feat_id = im_iter->i_;
			const IndexT second_feat_id = im_iter->j_;
			Vec2 first_feat = feats_1->GetRegionPosition(first_feat_id);
			Vec2 second_feat = feats_2->GetRegionPosition(second_feat_id);
			//check whether the first feature in first view is tracked.
			if (latter_opticaltable.count(Pair(first_view_id, first_feat_id)) == 0)
			{
				if (latter_view_id == first_view_id)
				{
					//the first feature  is not tracked in first view
					std::cerr << "Error: feature" << first_view_id << " is losed in view " << first_view_id << "\n";
					return false;
				}
				else
				{
					//the first feature  is not tracked in second view
					continue;
				}

			}
			//check whether the second feature in second view is tracked.
			if (latter_opticaltable.count(Pair(second_view_id, second_feat_id)) == 0)
			{
				if (latter_view_id == second_view_id)
				{
					//the second feature  is not tracked in second view
					std::cerr << "Error: feature" << second_feat_id << " is losed in view " << second_view_id << "\n";
					return false;
				}
				else
				{
					//the second feature  is not tracked in first view
					continue;
				}
			}
			//check the correctness of features
			const Optical_track& ot_1 = latter_opticaltable.at(Pair(first_view_id, first_feat_id));
			const Optical_track& ot_2 = latter_opticaltable.at(Pair(second_view_id, second_feat_id));

			if (double_comp(ot_1.rx, first_feat[0]) != 0 ||
				double_comp(ot_1.ry, first_feat[1]) != 0 ||
				double_comp(ot_2.rx, second_feat[0]) != 0 ||
				double_comp(ot_2.ry, second_feat[1]) != 0)
			{
				std::cerr << "Error: optical features is inconsistent with raw features\n";
				return false;
			}
			// compute the tracked distance of two features
			double dis_featpair = cal_dis(ot_1.x, ot_1.y, ot_2.x, ot_2.y);
			dis_featpairs.push_back(dis_featpair);

		}

	}

	//update distance threshold dynamically
	std::sort(dis_featpairs.begin(), dis_featpairs.end());
	int index_choose = (double)(dis_featpairs.size()) * 0.5;
	MaxDistanceThreshold_ = dis_featpairs[index_choose] * 10;
	return true;
}
  // the directory stores the binary file of features tracked in optical flow.
  std::string bin_dir_;
  // the distance threshold for optical filtering / matching
  double MaxDistanceThreshold_;
  // the repository storing position information  of features tracked by optical flow
  // map{view_id1 -> map{(view_id2,feat_id) -> track_infomation}}
  // meanning:
  //     *  the (`feat_id`)feature position on View(`view_id1`)  from View(`view_id2`) is recorded in track_infomation
  std::map<IndexT, std::map<Pair, Optical_track>> opticaltrack_table_;
  // the signal of successful reading
  bool readable_;

};



}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_CASCADE_HASHER_HPP
