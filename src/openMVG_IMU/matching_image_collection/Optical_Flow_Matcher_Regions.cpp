// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/matching_image_collection/Optical_Flow_Matcher_Regions.hpp"
#include "openMVG_IMU/utils/myoutput.hpp"
#include "openMVG_IMU/matching/cascade_hasher.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/metric.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/types.hpp"


#include "third_party/progress/progress.hpp"
#include <iomanip>
#include <fstream>
#include <iostream>

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

Optical_Flow_Matcher_Regions
::Optical_Flow_Matcher_Regions
(
	float distRatio, std::string bin_dir,double maxDistanceThreshold,const sfm::SfM_Data& sfm_data,
	std::string sMatches_dir
):Matcher(), f_dist_ratio_(distRatio), bin_dir_(bin_dir),MaxDistanceThreshold(maxDistanceThreshold),
sMatches_dir_(sMatches_dir)
{

	//read binary data of optical flow for every features
	std::cout << "**read binary features**\n";
	for (const auto& view_item : sfm_data.GetViews())
	{
		const IndexT view_id = view_item.first;
		std::stringstream ss;
		Optical_track ot;
		//ss << std::setw(5) << std::setfill('0') << std::stoi(stlplus::basename_part(view_item.second->s_Img_path)) << ".bin";
		ss << std::setw(5) << std::setfill('0') << view_id << ".bin";
		std::string filename = ss.str();
		std::ifstream file(stlplus::create_filespec(bin_dir, filename), std::ios::binary);

		if (!file.is_open())
		{
			std::cerr << stlplus::create_filespec(bin_dir, std::string(filename)) << " can not be read\n";
			return ;
		}
		if (opticaltrack_table.count(view_id))
		{
			std::cerr << "Error:duplicate view id\n";
			return ;
		}
		opticaltrack_table.emplace(view_id, std::map<Pair, Optical_track>());
		std::map<Pair, Optical_track>& map_optrack = opticaltrack_table.at(view_id);

		while (file.read(reinterpret_cast<char*>(&ot), sizeof(Optical_track)))
		{
			map_optrack[Pair(ot.view_id, ot.feat_id)] = ot;
			if (file.peek() == '\n')
			{
				file.ignore();
			}
		}
		
		std::cout << "#features of view  " << view_id << "("<< filename <<"): " << opticaltrack_table.at(view_id).size() << "\n";

	}

	
}

namespace impl
{

template <typename MatrixT, typename DistanceType>
bool Match_OpticalFlow
(
	const IndexT view_id_1,
	const HashedDescriptions& hashed_descriptions1,
	const MatrixT & descriptions1,
	const IndexT view_id_2,
	const HashedDescriptions& hashed_descriptions2,
	const MatrixT & descriptions2,
	const std::vector<features::PointFeature>& pts_1,
	const std::vector<features::PointFeature>& pts_2,
	const IndexT latter_view,
	const std::map<Pair, Optical_Flow_Matcher_Regions::Optical_track>& latter_opticaltable,
	const uint8_t nb_hash_code,
	IndMatches * pvec_indices,
	std::vector<DistanceType> * pvec_distances,
	const double MaxDistanceThreshold,
	const int NN = 2
)
{
	using MetricT = L2<typename MatrixT::Scalar>;
	MetricT metric;
	
	const int kNumTopCandidates = 10;
	const double eps = 1e-5;
	auto double_comp = [&](const double x, const double y)->bool {

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
	};

	auto cal_dis = [&](double x1, double y1, double x2, double y2)->double {
		return std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	};
		//optical flow match
		//Is the feature id consistent with index of vector<features::PointFeature> and column of descriptor matrix
		//

	
	Eigen::MatrixXi candidate_hamming_distances(
			hashed_descriptions2.hashed_desc.size(), nb_hash_code + 1);
	
	Eigen::VectorXi num_descriptors_with_hamming_distance(nb_hash_code + 1);
	std::vector<std::pair<DistanceType, int>> candidate_euclidean_distances;
	candidate_euclidean_distances.reserve(kNumTopCandidates);

	using HammingMetricType = matching::Hamming<stl::dynamic_bitset::BlockType>;
	static const HammingMetricType metricH = {};

	
	for (size_t i = 0; i < hashed_descriptions1.hashed_desc.size(); i++)
	{
		num_descriptors_with_hamming_distance.setZero();
		candidate_euclidean_distances.clear();
		size_t feats_selected = 0;
		// Compute the hamming distance of all feature pair whose distance is less than 10 px
		// Put the descriptors into buckets corresponding to their hamming distance.
		if (!latter_opticaltable.count(Pair(view_id_1, i)))
		{
			if (latter_view == view_id_1)
			{
				std::cout << "Error:feature " << i << "can not be found in the binary file of view " << view_id_1 << "\n";
				return false;
			}
			else
			{
				continue;
			}
		}
		const Optical_Flow_Matcher_Regions::Optical_track& ot_1 = latter_opticaltable.at(Pair(view_id_1, i));
		const auto& hashed_desc = hashed_descriptions1.hashed_desc[i];
		
		for (size_t j = 0; j < hashed_descriptions2.hashed_desc.size(); j++)
		{
			
			if (!latter_opticaltable.count(Pair(view_id_2, j)))
			{
				if (latter_view == view_id_2)
				{
					std::cout << "Error:feature " << j << "can not be found in the binary file of view " << view_id_2 << "\n";
					return false;
				}
				else
				{
					
					continue;
				}
			}
			
			const Optical_Flow_Matcher_Regions::Optical_track& ot_2 = latter_opticaltable.at(Pair(view_id_2, j));
			if (double_comp(ot_1.rx, pts_1[i].x()) != 0 ||
				double_comp(ot_1.ry, pts_1[i].y()) != 0 ||
				double_comp(ot_2.rx, pts_2[j].x()) != 0 ||
				double_comp(ot_2.ry, pts_2[j].y()) != 0)
			{
				std::cerr << "Error: optical features is inconsistent with raw features\n";
				return false;
			}
			
			double dis_featpair = cal_dis(ot_1.x, ot_1.y, ot_2.x, ot_2.y);
		
			if (dis_featpair <= MaxDistanceThreshold)
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


		// Compute the euclidean distance of the k descriptors with the best hamming
		// distance.
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
template <typename ScalarT>
bool Match
(
  const sfm::Regions_Provider & regions_provider,
  const Pair_Set & pairs,
  float fDistRatio,
  const std::string bin_dir,
  const std::string sMatches_dir,
  PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
  C_Progress * my_progress_bar,
	const double MaxDistanceThreshold,
	const std::map<IndexT, std::map<Pair, Optical_Flow_Matcher_Regions::Optical_track>>& opticaltrack_table
)
{
	
	////statistic////
	std::vector<size_t> vec_local_num_matches;
	std::vector<size_t> vec_num_matches_desfiltering;
	std::vector<size_t> vec_last_num_matches;


  if (!my_progress_bar)
    my_progress_bar = &C_Progress::dummy();
  my_progress_bar->restart(pairs.size(), "\n- Matching -\n");

  // Collect used view indexes
  std::set<IndexT> used_index;
  const int NN = 1;
  //Collect optical flow
  bool flag = true; //bc
  
  // Sort pairs according the first index to minimize later memory swapping
  using Map_vectorT = std::map<IndexT, std::vector<IndexT>>;
  Map_vectorT map_Pairs;
  for (const auto & pair_idx : pairs)
  {
    map_Pairs[pair_idx.first].push_back(pair_idx.second);
    used_index.insert(pair_idx.first);
    used_index.insert(pair_idx.second);
  }

  

  using BaseMat = Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  // Init the cascade hasher
  CascadeHasher_General cascade_hasher;
  if (!used_index.empty())
  {
    const IndexT I = *used_index.begin();
    const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
    const size_t dimension = regionsI->DescriptorLength();
    cascade_hasher.Init(dimension);
  }

  std::map<IndexT, HashedDescriptions> hashed_base_;

  // Compute the zero mean descriptor that will be used for hashing (one for all the image regions)
  Eigen::VectorXf zero_mean_descriptor;
  {
    Eigen::MatrixXf matForZeroMean;
    for (int i =0; i < used_index.size(); ++i)
    {
      std::set<IndexT>::const_iterator iter = used_index.begin();
      std::advance(iter, i);
      const IndexT I = *iter;
      const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
      const ScalarT * tabI =
        reinterpret_cast<const ScalarT*>(regionsI->DescriptorRawData());
      const size_t dimension = regionsI->DescriptorLength();
      if (i==0)
      {
        matForZeroMean.resize(used_index.size(), dimension);
        matForZeroMean.fill(0.0f);
      }
      if (regionsI->RegionCount() > 0)
      {
        Eigen::Map<BaseMat> mat_I( (ScalarT*)tabI, regionsI->RegionCount(), dimension);
        matForZeroMean.row(i) = CascadeHasher::GetZeroMeanDescriptor(mat_I);
      }
    }
    zero_mean_descriptor = CascadeHasher::GetZeroMeanDescriptor(matForZeroMean);
  }

  // Index the input regions
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i =0; i < used_index.size(); ++i)
  {
    std::set<IndexT>::const_iterator iter = used_index.begin();
    std::advance(iter, i);
    const IndexT I = *iter;
    const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
    const ScalarT * tabI =
      reinterpret_cast<const ScalarT*>(regionsI->DescriptorRawData());
    const size_t dimension = regionsI->DescriptorLength();

    Eigen::Map<BaseMat> mat_I( (ScalarT*)tabI, regionsI->RegionCount(), dimension);
#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
    {
      hashed_base_[I] =
        std::move(cascade_hasher.CreateHashedDescriptions(mat_I, zero_mean_descriptor));
    }
  }
  std::cout << "optical threshold: " << MaxDistanceThreshold << "\n";
  // Perform matching between all the pairs
  for (const auto & pair_it : map_Pairs)
  {
    if (my_progress_bar->hasBeenCanceled())
      break;
    const IndexT I = pair_it.first;
    const std::vector<IndexT> & indexToCompare = pair_it.second;

    const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
    if (regionsI->RegionCount() == 0)
    {
      (*my_progress_bar) += indexToCompare.size();
      continue;
    }

    const std::vector<features::PointFeature> pointFeaturesI = regionsI->GetRegionsPositions();
    const ScalarT * tabI =
      reinterpret_cast<const ScalarT*>(regionsI->DescriptorRawData());
    const size_t dimension = regionsI->DescriptorLength();
    Eigen::Map<BaseMat> mat_I( (ScalarT*)tabI, regionsI->RegionCount(), dimension);

#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      if (my_progress_bar->hasBeenCanceled())
        continue;
	  
      const size_t J = indexToCompare[j];
      const std::shared_ptr<features::Regions> regionsJ = regions_provider.get(J);
	
      if (regionsI->Type_id() != regionsJ->Type_id())
      {
        ++(*my_progress_bar);
        continue;
      }
	  ////statistic////
	  size_t local_num_matches = 0;
	  size_t num_matches_desfiltering = 0;
	  size_t last_num_matches = 0;

      // Matrix representation of the query input data;
      const ScalarT * tabJ = reinterpret_cast<const ScalarT*>(regionsJ->DescriptorRawData());
      Eigen::Map<BaseMat> mat_J( (ScalarT*)tabJ, regionsJ->RegionCount(), dimension);

	  ////find features in global
      IndMatches pvec_indices_global;
      using ResultType = typename Accumulator<ScalarT>::Type;
      std::vector<ResultType> pvec_distances_global;
	  pvec_distances_global.reserve(regionsJ->RegionCount() * 2);
	  pvec_indices_global.reserve(regionsJ->RegionCount() * 2);
	  const std::vector<features::PointFeature> pointFeaturesJ = regionsJ->GetRegionsPositions();
	  cascade_hasher.Match_HashedDescriptions<BaseMat, ResultType>(
		  hashed_base_[J], mat_J,
		  hashed_base_[I], mat_I,
		  &pvec_indices_global, &pvec_distances_global,1);


	  
	  //match by optical flow
	  //std::cout << "Pair " << I << "-" << J << "\n";

	  //find features in local
	  IndMatches pvec_indices_local;
	  std::vector<ResultType> pvec_distances_local;
	  pvec_distances_local.reserve(regionsJ->RegionCount() * 2);
	  pvec_indices_local.reserve(regionsJ->RegionCount() * 2);
	  if(!Match_OpticalFlow<BaseMat, ResultType>
		  (
			  J,
			  hashed_base_[J], mat_J,
			  I,
			  hashed_base_[I], mat_I,

			  pointFeaturesJ,
			  pointFeaturesI,
			  (I >= J ? I : J),
			  opticaltrack_table.at((I >= J ? I : J)),
			  cascade_hasher.nb_hash_code_,               //cascade_hasher is initialized with the dimension
			  &pvec_indices_local,
			  &pvec_distances_local,
			  MaxDistanceThreshold,
			  1
			  ))
	  {
		  std::cout << "Error:in matching\n";
		  assert(0);
	  }
	  local_num_matches = pvec_indices_local.size();
	  //
	  //std::cout << "number of candidate feat pair found:" << pvec_indices.size()/2 << "\n";
	  
      // Match the query descriptors to the database
      /*cascade_hasher.Match_HashedDescriptions<BaseMat, ResultType>(
        hashed_base_[J], mat_J,
        hashed_base_[I], mat_I,
        &pvec_indices, &pvec_distances);*/
	  
	  matching::IndMatches vec_putative_matches;
	  vec_putative_matches.reserve(pvec_indices_local.size());
	  {
		  const float MinDesDistance = 1.1;

		  for (size_t ki = 0; ki < pvec_indices_local.size(); ki++)
		  {
			  for (size_t wi = 0; wi < pvec_indices_global.size(); wi++)
			  {
				  if (pvec_indices_local[ki].i_ == pvec_indices_global[wi].i_ &&
					  pvec_distances_local[ki] < pvec_distances_global[wi] * MinDesDistance)
				  {
					  vec_putative_matches.emplace_back(pvec_indices_local[ki].j_, pvec_indices_local[ki].i_);
				  }
			  }
		  }
	  }
	  num_matches_desfiltering = vec_putative_matches.size();
	  ///openmvg s///
		  // // Filter the matches using a distance ratio test:
		  // //   The probability that a match is correct is determined by taking
		  // //   the ratio of distance from the closest neighbor to the distance
		  // //   of the second closest.
		  // matching::NNdistanceRatio(
			 //pvec_distances.begin(), // distance start
			 //pvec_distances.end(),   // distance end
			 //NN, // Number of neighbor in iterator sequence (minimum required 2)
			 //vec_nn_ratio_idx, // output (indices that respect the distance Ratio)
			 //  1.1);
		  // 
		  // 
		  // for (size_t k=0; k < vec_nn_ratio_idx.size(); ++k)
		  // {
			 //const size_t index = vec_nn_ratio_idx[k];
			 //vec_putative_matches.emplace_back(pvec_indices[index*2].j_, pvec_indices[index*2].i_);
		  // }


	 
	  ///openmvg e///
	 // ////bc s///

	 // matching::IndMatches vec_putative_matches;
		//vec_putative_matches.reserve(pvec_indices.size());
		//for (size_t k=0; k < pvec_indices.size(); ++k)
		//{
		//	vec_putative_matches.emplace_back(pvec_indices[k].j_, pvec_indices[k].i_);
		//}

	 // ////bc e///
      // Remove duplicates
      matching::IndMatch::getDeduplicated(vec_putative_matches);

      // Remove matches that have the same (X,Y) coordinates
	  
      matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
        pointFeaturesI, pointFeaturesJ);
      matchDeduplicator.getDeduplicated(vec_putative_matches);
	  
	  last_num_matches = vec_putative_matches.size();
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
      {
		  std::cout << "Pair " << I << "-" << J << "\n";
		  std::cout << "number of candidate feat pair found:" << vec_putative_matches.size() << "\n";

		  vec_local_num_matches.emplace_back(local_num_matches);
		  vec_num_matches_desfiltering.emplace_back(num_matches_desfiltering);
		  vec_last_num_matches.emplace_back(last_num_matches);
        if (!vec_putative_matches.empty())
        {
          map_PutativesMatches.insert(
            {
              {I,J},
              std::move(vec_putative_matches)
            });
        }
      }
      ++(*my_progress_bar);
    }
  }
  {
	  std::cout << "The statistic of local_num_matches:\n";
	  sort(vec_local_num_matches.begin(), vec_local_num_matches.end());
	  utils::MinMaxMedianMean<size_t>(std::cout, vec_local_num_matches);
	  std::ofstream local_num_matches_log(stlplus::create_filespec(sMatches_dir, "local_num_matches.txt"));
	  utils::MinMaxMedianMean<size_t>(local_num_matches_log, vec_local_num_matches);
	  local_num_matches_log.close();

	  std::cout << "The statistic of desfiltering_num_matches:\n";
	  sort(vec_num_matches_desfiltering.begin(), vec_num_matches_desfiltering.end());
	  utils::MinMaxMedianMean<size_t>(std::cout, vec_num_matches_desfiltering);
	  std::ofstream desfiltering_num_matches_log(stlplus::create_filespec(sMatches_dir, "desfiltering_num_matches.txt"));
	  utils::MinMaxMedianMean<size_t>(desfiltering_num_matches_log, vec_num_matches_desfiltering);
	  desfiltering_num_matches_log.close();

	  std::cout << "The statistic of last_num_matches:\n";
	  sort(vec_last_num_matches.begin(), vec_last_num_matches.end());
	  utils::MinMaxMedianMean<size_t>(std::cout, vec_last_num_matches);
	  std::ofstream last_num_matches_log(stlplus::create_filespec(sMatches_dir, "last_num_matches.txt"));
	  utils::MinMaxMedianMean<size_t>(last_num_matches_log, vec_last_num_matches);
	  last_num_matches_log.close();
  }
  return flag;
}



} // namespace impl

void Optical_Flow_Matcher_Regions::Match
(
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
  C_Progress * my_progress_bar
)const
{
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
  if (!regions_provider)
    return;

  if (regions_provider->IsBinary())
    return;

  if (regions_provider->Type_id() == typeid(unsigned char).name())
  {
    impl::Match<unsigned char>(
      *regions_provider.get(),
      pairs,
      f_dist_ratio_,bin_dir_,
		sMatches_dir_,
      map_PutativesMatches,
      my_progress_bar,
		MaxDistanceThreshold,
		opticaltrack_table);
  }
  else
  if (regions_provider->Type_id() == typeid(float).name())
  {
	  impl::Match<float>(
      *regions_provider.get(),
      pairs,
      f_dist_ratio_, bin_dir_,
		  sMatches_dir_,
      map_PutativesMatches,
      my_progress_bar,
		  MaxDistanceThreshold,
		  opticaltrack_table);
  }
  else
  {
    std::cerr << "Matcher not implemented for this region type" << std::endl;
  }
}



} // namespace openMVG
} // namespace matching_image_collection
