// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/matching_image_collection/Hierarchical_Matcher_Regions.hpp"

#include "openMVG_IMU/matching/cascade_hasher2.hpp"
#include "openMVG_IMU/utils/myoutput.hpp"
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
#include <typeinfo>

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

Hierarchical_Matcher_Regions
::Hierarchical_Matcher_Regions
(
	float distRatio, std::string bin_dir, double maxDistanceThreshold, const sfm::SfM_Data& sfm_data,
	const std::string sMatches_dir, bool bfeature_validation, bool bopticalfiltering,
    bool bdynamicdistance, bool bopticalmatching,bool bdebug
):Matcher(), f_dist_ratio_(distRatio),sMatches_dir_(sMatches_dir),bdebug_(bdebug),
opticalflow_container(bin_dir,maxDistanceThreshold,sfm_data),
bfeature_validation_(bfeature_validation), bopticalfiltering_(bopticalfiltering),
bdynamicdistance_(bdynamicdistance), bopticalmatching_(bopticalmatching)
{
	
	
}

namespace impl
{


template <typename ScalarT>
void Match_Hierarchical
(
	const sfm::Regions_Provider & regions_provider,
	const Pair_Set & pairs,
	float fDistRatio,
	PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
	OpticalFlow_Container opticalflow_container,
	const std::string sMatchesDir,
	C_Progress * my_progress_bar,
	std::map<IndexT,std::set<IndexT>>& notable_features,
	bool bfeature_validation = true,
	bool bopticalfiltering = true,
	bool bdynamicdistance = true,
	bool bopticalmatching = true,
	bool bnotablevalidation = false,         
	bool bdebug = true
)
{
  
  using ResultType = typename Accumulator<ScalarT>::Type;
  ////statistic////
  std::vector<size_t> first_num_matches;
  std::vector<size_t> second_num_matches;
  std::vector<size_t> third_num_matches;
  std::vector<size_t> fourth_num_matches;
  std::vector<size_t> fifth_num_matches;
  std::vector<size_t> sixth_num_matches;
  std::vector<size_t> seventh_num_matches;
  std::vector<size_t> eighth_num_matches;

  std::map<Pair, std::map<Pair, ResultType>> descdis_putativematches;
  ResultType MaxDescDistance;
  if (typeid(MaxDescDistance).name() == typeid(float).name())
  {
	  std::cout << "ResultType:float\n";
	  MaxDescDistance = 0.0;
  }
  else
  {
	  std::cout << "Error:Unknown ResultType(" << typeid(MaxDescDistance).name() << "\n";
	  return;
  }
  if (!my_progress_bar)
	  my_progress_bar = &C_Progress::dummy();
  my_progress_bar->restart(pairs.size(), "\n- Matching -\n");
  ///////////////
  // Collect used view indexes
  std::set<IndexT> used_index;
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
  CascadeHasher2 cascade_hasher2;
  if (!used_index.empty())
  {
    const IndexT I = *used_index.begin();
    const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
    const size_t dimension = regionsI->DescriptorLength();
    cascade_hasher2.Init(dimension);
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
        std::move(cascade_hasher2.CreateHashedDescriptions(mat_I, zero_mean_descriptor));
    }
  }
  
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
	  size_t num_matches_first = 0;
	  size_t num_matches_second = 0;

      // Matrix representation of the query input data;
      const ScalarT * tabJ = reinterpret_cast<const ScalarT*>(regionsJ->DescriptorRawData());
      Eigen::Map<BaseMat> mat_J( (ScalarT*)tabJ, regionsJ->RegionCount(), dimension);

      IndMatches pvec_indices_JI;
      
      std::vector<ResultType> pvec_distances_JI;
      pvec_distances_JI.reserve(regionsJ->RegionCount() * 2);
      pvec_indices_JI.reserve(regionsJ->RegionCount() * 2);
	  //find one feature in I for every feature in J
      // Match the query descriptors to the database
      cascade_hasher2.Match_HashedDescriptions<BaseMat, ResultType>(
        hashed_base_[J], mat_J,
        hashed_base_[I], mat_I,
        &pvec_indices_JI, &pvec_distances_JI);

      std::vector<int> vec_nn_ratio_idx_JI;
      // Filter the matches using a distance ratio test:
      //   The probability that a match is correct is determined by taking
      //   the ratio of distance from the closest neighbor to the distance
      //   of the second closest.
      matching::NNdistanceRatio(
        pvec_distances_JI.begin(), // distance start
        pvec_distances_JI.end(),   // distance end
        2, // Number of neighbor in iterator sequence (minimum required 2)
        vec_nn_ratio_idx_JI, // output (indices that respect the distance Ratio)
        Square(fDistRatio));
	 num_matches_first = vec_nn_ratio_idx_JI.size();
	IndMatches pvec_indices_IJ;
	std::vector<int> vec_nn_ratio_idx_IJ;
	
	//std::cout << "#first mathcing(" << J << "->" << I << "): " << vec_nn_ratio_idx_JI.size() << " feature pairs are found.\n";
	//find one feature in J for every feature in I
	if(bfeature_validation)
	{
		std::set<IndexT> matched_featIids;
		std::vector<ResultType> pvec_distances_IJ;
		
		for (size_t k=0; k < vec_nn_ratio_idx_JI.size(); ++k)
		{
			const size_t index = vec_nn_ratio_idx_JI[k];
			matched_featIids.insert(pvec_indices_JI[index*2].j_);
		}
		pvec_distances_IJ.reserve(matched_featIids.size() * 2);
		pvec_indices_IJ.reserve(matched_featIids.size() * 2);
		cascade_hasher2.Match_HashedDescriptions_Specifiedfeat1<BaseMat, ResultType>(
		hashed_base_[I], mat_I,
		hashed_base_[J], mat_J,
		&pvec_indices_IJ, &pvec_distances_IJ,
		matched_featIids, 2);

			matching::NNdistanceRatio(
		pvec_distances_IJ.begin(), // distance start
		pvec_distances_IJ.end(),   // distance end
		2, // Number of neighbor in iterator sequence (minimum required 2)
		vec_nn_ratio_idx_IJ, // output (indices that respect the distance Ratio)
		Square(fDistRatio));

	}
	const std::vector<features::PointFeature> pointFeaturesJ = regionsJ->GetRegionsPositions();

	//check the feautre pair both appear in two matches.
	matching::IndMatches vec_putative_matches;
	std::map<Pair, ResultType> map_putative_descdis;
	bool bindex_consistent = true;
	vec_putative_matches.reserve(vec_nn_ratio_idx_JI.size());
	for (size_t k = 0; k < vec_nn_ratio_idx_JI.size(); ++k)
	{
		const size_t index1 = vec_nn_ratio_idx_JI[k];
		const IndexT feat_i_id1 = pvec_indices_JI[index1 * 2].j_;
		const IndexT feat_j_id1 = pvec_indices_JI[index1 * 2].i_;
		const features::PointFeature pointfeati_1 = pointFeaturesI[feat_i_id1];
		const features::PointFeature pointfeatj_1 = pointFeaturesJ[feat_j_id1];
		bool boccur_featurepair = false;
		if (!bfeature_validation)
		{
			boccur_featurepair = true;
		}
		else
		{

			for (size_t w = 0; w < vec_nn_ratio_idx_IJ.size(); ++w)
			{
				const size_t index2 = vec_nn_ratio_idx_IJ[w];
				const IndexT feat_i_id2 = pvec_indices_IJ[index2 * 2].i_;
				const features::PointFeature pointfeati_2 = pointFeaturesI[feat_i_id2];
				if (feat_i_id2 != feat_i_id1 && pointfeati_2.coords() != pointfeati_1.coords()) continue;

				const IndexT feat_j_id2 = pvec_indices_IJ[index2 * 2].j_;
				
				const features::PointFeature pointfeatj_2 = pointFeaturesJ[feat_j_id2];
				if (pointfeati_2.coords() != pointfeati_1.coords())
				{
					std::cout << "Error:index of features in view I is not aligned\n";
					bindex_consistent = false;
					break;
				}
				if (feat_j_id1 == feat_j_id2 || (pointfeatj_1.coords() == pointfeatj_2.coords()))
				{
					boccur_featurepair = true;
					break;
				}
				//break;
			}
		}
		if (!bindex_consistent) break;
		if (boccur_featurepair)
		{
			vec_putative_matches.emplace_back(feat_i_id1, feat_j_id1);
			map_putative_descdis.emplace(Pair(feat_i_id1, feat_j_id1), pvec_distances_JI[index1 * 2]);
		}

	}
	
	num_matches_second = vec_putative_matches.size();
	//std::cout<<"#mathcing("<<I<<"->"<<J<<"): "<<vec_putative_matches.size()<<" feature pairs are found.\n";


  // Remove duplicates
	matching::IndMatch::getDeduplicated(vec_putative_matches);

	// Remove matches that have the same (X,Y) coordinates

	matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
		pointFeaturesI, pointFeaturesJ);
	matchDeduplicator.getDeduplicated(vec_putative_matches);



#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
	{
		//assert(bindex_consistent);   //?
		if (!bindex_consistent)
		{
			std::cout << "View " << I << "-" << J << " :index consistent\n";

		}
		if (!vec_putative_matches.empty())
		{

			map_PutativesMatches.insert(
			{
			  {I,J},
			  std::move(vec_putative_matches)
			});
			descdis_putativematches.emplace(Pair(I, J), std::move(map_putative_descdis));
			first_num_matches.push_back(num_matches_first);
			second_num_matches.push_back(num_matches_second);

		}
	}
	++(*my_progress_bar);
	}
  }
  std::cout << "The statistic of first_num_matches:\n";
  sort(first_num_matches.begin(), first_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, first_num_matches);
  std::cout << "The statistic of second_num_matches:\n";
  sort(second_num_matches.begin(), second_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, second_num_matches);
  if (!Save(map_PutativesMatches,
	  std::string(sMatchesDir + "/" + "matches.siftfeat.bin")))
  {
	  std::cerr
		  << "Cannot save computed matches in: "
		  << std::string(sMatchesDir + "/" + "matches.siftfeat.bin");
  }
  if (bdebug)
  {
	  std::ofstream third_num_matches_log(stlplus::create_filespec(sMatchesDir, "third_num_matches.txt"));
	  utils::Output_Matches(third_num_matches_log, map_PutativesMatches, true);
	  third_num_matches_log.close();
  }
	  //////////////optical filtering/////////////////
	if (!bopticalfiltering)
	{
		return;
	}
	if (!opticalflow_container.readable_)
	{
		std::cout << "the binary files are not read completely.\n";
		return;
	}
	if (bdynamicdistance)
	{
		bool succ = opticalflow_container.SetDynamicDistanceThreshold(map_PutativesMatches, &regions_provider);
		if (succ)
		{
			std::cout << "dynamic threshold is set successfully as " << opticalflow_container.MaxDistanceThreshold_ << "\n";
		}
		else
		{
			std::cout << "dynamic threshold is set failed,use static threshold:" << opticalflow_container.MaxDistanceThreshold_ << "\n";
		}
	}
	  
	  
	  
	  std::vector<double> vec_ofdis;
	 bool bfiltering = opticalflow_container.Optical_Filtering
							(map_PutativesMatches, &regions_provider,my_progress_bar,
							third_num_matches,fourth_num_matches, vec_ofdis
							);

	 if (!bfiltering)
	 {
		 return;
	  }
	 {
		 /////compute notable features in every frame
		 //map{{view_id}->map{{feature_id}->{status}}}
		 std::map<int, std::map<int, int>>  notablefeature_records;
		 const IndexT NumAdjFrames = 2;
		 for (const auto& pairwise_item : map_PutativesMatches)
		 {
			 const IndexT I = pairwise_item.first.first;
			 const IndexT J = pairwise_item.first.second;
			 if (J - I > NumAdjFrames) continue;   //3 adjacent frames
			 if (!notablefeature_records.count(I))
			 {
				 notablefeature_records.emplace(I, std::map<int, int>());
			 }

			 for (const auto& indmatch : pairwise_item.second)
			 {
				 if (!notablefeature_records.at(I).count(indmatch.i_))
				 {
					 notablefeature_records.at(I).emplace(indmatch.i_, 0);
				 }
				 int& mark = notablefeature_records.at(I).at(indmatch.i_);
				 mark |= 1 << (J - I - 1);
			 }
		 }

		 for (const auto& features_in_view : notablefeature_records)
		 {
			 notable_features.emplace(features_in_view.first, std::set<IndexT>());
			 for (const auto& feature_items : features_in_view.second)
			 {
				 if (feature_items.second == 3)
					 notable_features.at(features_in_view.first).insert(feature_items.first);
			 }
		 }
		 if (bnotablevalidation)
		 {
			 std::cout << "notable validation for 3 adjacent frames\n";
			 PairWiseMatches::iterator pm_iter;
			 size_t num_deleted = 0;
			 for (pm_iter = map_PutativesMatches.begin(); pm_iter != map_PutativesMatches.end(); pm_iter++)
			 {
				 const IndexT I = pm_iter->first.first;
				 const IndexT J = pm_iter->first.second;
				 if (J - I > NumAdjFrames) continue;   //3 adjacent frames
				 IndMatches::iterator im_iter;
				 if (!notable_features.count(I) || !notable_features.count(J))
				 {
					 continue;
				 }
				 for (im_iter = pm_iter->second.begin(); im_iter != pm_iter->second.end();)
				 {
					 if (notable_features.at(I).count(im_iter->i_) && notable_features.at(J).count(im_iter->j_))
					 {
						 im_iter++;
					 }
					 else
					 {
						 num_deleted++;
						 im_iter = pm_iter->second.erase(im_iter);
					 }
				 }
			 }
			 std::cout << num_deleted << " are deleted.\n";
		 }
		 
	 }
	  {
		  if (bdebug) {
			  std::ofstream ofdis_log(stlplus::create_filespec(sMatchesDir, "ofdis_log.txt"));
			  std::sort(vec_ofdis.begin(), vec_ofdis.end());
			  utils::AllMinMaxMedianMean<double>(ofdis_log, vec_ofdis);
			  ofdis_log.close();
		  }
		  
		  //delete the distance records of filtered feature pair
		  std::vector<ResultType> vec_descdis;
		  std::ofstream descdis_pairlog(stlplus::create_filespec(sMatchesDir, "descdis_pairlog.txt"));
		  for (const auto& pairwise_item : map_PutativesMatches)
		  {

			  std::set<Pair> existed_featpair;
			  size_t cur_size = vec_descdis.size();
			  for (const auto& indmatch : pairwise_item.second)
			  {
				  existed_featpair.insert(Pair(indmatch.i_, indmatch.j_));
			  }
			  descdis_pairlog << "View Pair " << pairwise_item.first.first << "," << pairwise_item.first.second << "\n";
			  std::map<Pair, ResultType>& map_descdis_IJ = descdis_putativematches.at(pairwise_item.first);
			  for (std::map<Pair, ResultType>::iterator it = map_descdis_IJ.begin(); it != map_descdis_IJ.end();)
			  {
				  if (!existed_featpair.count(it->first))
				  { 
					  it = map_descdis_IJ.erase(it);
				  }
				  else
				  {
					  MaxDescDistance = MaxDescDistance > it->second ? MaxDescDistance : it->second;
					  vec_descdis.emplace_back(it->second);
					  //descdis_pairlog << "	" << it->first.first << "-" << it->first.second <<":"<< it->second <<"\n";
					  it++;
				  }
			  }
			  std::vector<ResultType> vec_descdis_local(vec_descdis.begin() + cur_size, vec_descdis.end());
			  std::sort(vec_descdis_local.begin(), vec_descdis_local.end());
			  utils::MinMaxMedianMean(descdis_pairlog, vec_descdis_local);

		  }
		  if (bdebug)
		  {
			  std::ofstream descdis_log(stlplus::create_filespec(sMatchesDir, "desdis_log.txt"));
			  std::sort(vec_descdis.begin(), vec_descdis.end());
			  utils::AllMinMaxMedianMean<ResultType>(descdis_log, vec_descdis);
			  std::ofstream fourth_num_matches_log(stlplus::create_filespec(sMatchesDir, "fourth_num_matches.txt"));
			  utils::Output_Matches(fourth_num_matches_log, map_PutativesMatches, true);
			  descdis_log.close();
			  fourth_num_matches_log.close();
		  }
		  
		 
		  descdis_pairlog.close();
		  
		  
		  std::cout << "The statistic of third_num_matches:\n";
		  sort(third_num_matches.begin(), third_num_matches.end());
		  utils::MinMaxMedianMean<size_t>(std::cout, third_num_matches);
		  std::cout << "The statistic of fourth_num_matches:\n";
		  sort(fourth_num_matches.begin(), fourth_num_matches.end());
		  utils::MinMaxMedianMean<size_t>(std::cout, fourth_num_matches);
	  }
	  //////////////optical matching/////////////////
	  if (!bopticalmatching)
	  {
		  return;
	  }
	  std::cout << "optical threshold: " << opticalflow_container.MaxDistanceThreshold_ << "\n";
	  typedef std::pair<IndexT, std::vector<IndMatch>*> Pairiva;
	  std::map<IndexT, std::vector<Pairiva>> map_putmatches;

	  for (auto& putmatch_item : map_PutativesMatches)
	  {
		  map_putmatches[putmatch_item.first.first].emplace_back(
			  Pairiva(putmatch_item.first.second, &putmatch_item.second));
	  }
	  if (!my_progress_bar)
		  my_progress_bar = &C_Progress::dummy();
	  my_progress_bar->restart(map_PutativesMatches.size(), "\n- Optical Matching -\n");
	  //MaxDescDistance = MaxDescDistance;
	  std::cout<<"dynamic maximum descriptor :"<<MaxDescDistance<<"\n";
	  for(auto& putmatch_item: map_putmatches)
	  {
		  if (my_progress_bar->hasBeenCanceled())
			  break;
		  	const IndexT I = putmatch_item.first;
			
			const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
			const std::vector<features::PointFeature> pointFeaturesI = regionsI->GetRegionsPositions();
			if (regionsI->RegionCount() == 0)
			{
				(*my_progress_bar) += putmatch_item.second.size();
				continue;
			}
			const ScalarT * tabI =
				reinterpret_cast<const ScalarT*>(regionsI->DescriptorRawData());
			const size_t dimension = regionsI->DescriptorLength();
			Eigen::Map<BaseMat> mat_I((ScalarT*)tabI, regionsI->RegionCount(), dimension);
			

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
			for(int j = 0 ; j< (int)putmatch_item.second.size(); j++)
			{
				if (my_progress_bar->hasBeenCanceled())
					continue;

				
				
				Pairiva& pairiva = putmatch_item.second[j];
				const IndexT J = pairiva.first;
				//find dynamic similarity threshold for image pair
				ResultType MaxSimilarity_perpair = 0.0;
				for (const auto& descdis_item : descdis_putativematches.at(Pair(I, J)))
				{
					MaxSimilarity_perpair = MaxSimilarity_perpair > descdis_item.second ?
						MaxSimilarity_perpair : descdis_item.second;
				}

				///static////
				size_t num_fourth_matches = pairiva.second->size();
				size_t num_fifth_matches = 0;
				size_t num_sixth_matches = 0;
				
				std::set<IndexT> paired_featiid;
				std::set<IndexT> paired_featjid;
				const std::shared_ptr<features::Regions> regionsJ = regions_provider.get(J);
				const std::vector<features::PointFeature> pointFeaturesJ = regionsJ->GetRegionsPositions();

				const ScalarT * tabJ = reinterpret_cast<const ScalarT*>(regionsJ->DescriptorRawData());
				
				Eigen::Map<BaseMat> mat_J((ScalarT*)tabJ, regionsJ->RegionCount(), dimension);
				
				
				if (regionsI->Type_id() != regionsJ->Type_id())
				{
					++(*my_progress_bar);
					continue;
				}
				////find features in global
				IndMatches pvec_indices_global;
				
				std::vector<ResultType> pvec_distances_global;
				pvec_distances_global.reserve(regionsJ->RegionCount() * 2);
				pvec_indices_global.reserve(regionsJ->RegionCount() * 2);
				
				cascade_hasher2.Match_HashedDescriptions<BaseMat, ResultType>(
					hashed_base_[J], mat_J,
					hashed_base_[I], mat_I,
					&pvec_indices_global, &pvec_distances_global, 1);

				for (const auto& featmatch_item : (*pairiva.second))
				{
					paired_featiid.insert(featmatch_item.i_);
					paired_featjid.insert(featmatch_item.j_);
				}
				std::vector<ResultType> pvec_distances_JI;
				IndMatches pvec_indices_JI;
				bool bofmatching=opticalflow_container.Optical_Matching<BaseMat, ResultType>
					(
						J,
						hashed_base_[J], mat_J,
						I,
						hashed_base_[I], mat_I,

						pointFeaturesJ,
						pointFeaturesI,
						paired_featjid,
						paired_featiid,
						(I >= J ? I : J),

						cascade_hasher2.nb_hash_code_,               //cascade_hasher is initialized with the dimension
						&pvec_indices_JI,
						&pvec_distances_JI,
						1
						);
				num_fifth_matches = pvec_indices_JI.size() + num_fourth_matches;
				matching::IndMatches vec_putative_matches;
				// set threshold  as MaxSimilarity_perpair
				/*for (size_t k = 0; k < pvec_indices_JI.size(); ++k)
				{
					if (pvec_distances_JI[k] < MaxSimilarity_perpair)
						vec_putative_matches.emplace_back(pvec_indices_JI[k].j_, pvec_indices_JI[k].i_);
				}*/
				vec_putative_matches.reserve(pvec_indices_JI.size());
				{
					const float MinDesDistance = 1.1;

					for (size_t ki = 0; ki < pvec_indices_JI.size(); ki++)
					{
						for (size_t wi = 0; wi < pvec_indices_global.size(); wi++)
						{
							if (pvec_indices_JI[ki].i_ == pvec_indices_global[wi].i_ &&
								pvec_distances_JI[ki] < pvec_distances_global[wi] * MinDesDistance)
							{
								vec_putative_matches.emplace_back(pvec_indices_JI[ki].j_, pvec_indices_JI[ki].i_);
							}
						}
					}
				}
				num_sixth_matches = vec_putative_matches.size() + num_fourth_matches;
				// Remove duplicates
				matching::IndMatch::getDeduplicated(vec_putative_matches);

				// Remove matches that have the same (X,Y) coordinates

				matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
					pointFeaturesI, pointFeaturesJ);
				matchDeduplicator.getDeduplicated(vec_putative_matches);
				
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif		
				{
					//assert(bofmatching);
					/*std::cout << "fifth num(:" << I << "-" << J << "):" << num_fifth_matches - num_fourth_matches << "\n";
					
					std::cout << "descriptor filtering num(:" << I << "-" << J << "):" << num_sixth_matches - num_fourth_matches << "\n";
					
					std::cout << "seven num(:" << I << "-" << J << "):" << vec_putative_matches.size() << "\n";*/
					
					pairiva.second->insert(pairiva.second->end(), vec_putative_matches.begin(), vec_putative_matches.end());
					/*std::cout << "fourth num all(:" << I << "-" << J << "):" << num_fourth_matches << "\n";
					std::cout << "sixth num all(:" << I << "-" << J << "):" << num_sixth_matches << "\n";
					std::cout << "fifth num all(:" << I << "-" << J << "):" << num_fifth_matches << "\n";
					std::cout << "*final num all(:" << I << "-" << J << "):" << pairiva.second->size() << "\n";*/
					fifth_num_matches.emplace_back(num_fifth_matches);
					sixth_num_matches.emplace_back(num_sixth_matches);
					seventh_num_matches.emplace_back(pairiva.second->size());
					// Remove duplicates
					matching::IndMatch::getDeduplicated(*pairiva.second);

					// Remove matches that have the same (X,Y) coordinates

					matching::IndMatchDecorator<float> matchDeduplicator(*pairiva.second,
						pointFeaturesI, pointFeaturesJ);
					matchDeduplicator.getDeduplicated(*pairiva.second);
					eighth_num_matches.emplace_back(pairiva.second->size());
				}
				++(*my_progress_bar);

			}
			
			
	  }
	  
	
  std::cout<<"The statistic of fifth_num_matches:\n";
  sort(fifth_num_matches.begin(), fifth_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout,fifth_num_matches);
  std::cout << "The statistic of sixth_num_matches:\n";
  sort(sixth_num_matches.begin(), sixth_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, sixth_num_matches);
  std::cout << "The statistic of seventh_num_matches:\n";
  sort(seventh_num_matches.begin(), seventh_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, seventh_num_matches);
  std::cout << "The statistic of eighth_num_matches:\n";
  sort(eighth_num_matches.begin(), eighth_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, eighth_num_matches);
  if(bdebug)
  {
	  std::ofstream fifth_num_matches_log(stlplus::create_filespec(sMatchesDir, "final_matches.txt"));
	  utils::Output_Matches(fifth_num_matches_log, map_PutativesMatches,true);
  }
  
  
}
} // namespace impl



void Hierarchical_Matcher_Regions::Hierarchical_Match
(const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
	const Pair_Set & pairs,
	matching::PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
	C_Progress * my_progress_bar
)
{
#ifdef OPENMVG_USE_OPENMP
	std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
	if (!regions_provider)
		return;

	if (regions_provider->IsBinary())
		return;
	/*bool bfeature_validation = false;
	bool bopticalfiltering = false;
	bool bdynamicdistance = false;
	bool bopticalmatching = false;*/
	bool bnotablevalidation = true;   //turn it down
	if (regions_provider->Type_id() == typeid(unsigned char).name())
	{
		impl::Match_Hierarchical<unsigned char>(
			*regions_provider.get(),
			pairs,
			f_dist_ratio_,
			map_PutativesMatches,
			opticalflow_container,
			sMatches_dir_,
			my_progress_bar, 
			notable_features,
			bfeature_validation_,
			bopticalfiltering_,
			bdynamicdistance_,
			bopticalmatching_, bnotablevalidation,bdebug_);
	}
	else if (regions_provider->Type_id() == typeid(float).name())
	{
		impl::Match_Hierarchical<float>(
			*regions_provider.get(),
			pairs,
			f_dist_ratio_,
			map_PutativesMatches,
			opticalflow_container,
			sMatches_dir_,
			my_progress_bar,
			notable_features, 
			bfeature_validation_,
			bopticalfiltering_,
			bdynamicdistance_,
			bopticalmatching_, bnotablevalidation,bdebug_);
	}
	else
	{
		std::cerr << "Matcher not implemented for this region type" << std::endl;
		return;
	}
	

}

void Hierarchical_Matcher_Regions::Match
(
	const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
	const Pair_Set & pairs,
	PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
	C_Progress * my_progress_bar
)const
{
//#ifdef OPENMVG_USE_OPENMP
//	std::cout << "Using the OPENMP thread interface" << std::endl;
//#endif
//	if (!regions_provider)
//		return;
//
//	if (regions_provider->IsBinary())
//		return;
//
//	if (regions_provider->Type_id() == typeid(unsigned char).name())
//	{
//		impl::Match<unsigned char>(
//			*regions_provider.get(),
//			pairs,
//			f_dist_ratio_, bin_dir_,
//			map_PutativesMatches,
//			my_progress_bar,
//			MaxDistanceThreshold,
//			opticaltrack_table);
//	}
//	else
//		if (regions_provider->Type_id() == typeid(float).name())
//		{
//			impl::Match<float>(
//				*regions_provider.get(),
//				pairs,
//				f_dist_ratio_, bin_dir_,
//				map_PutativesMatches,
//				my_progress_bar,
//				MaxDistanceThreshold,
//				opticaltrack_table);
//		}
//		else
//		{
//			std::cerr << "Matcher not implemented for this region type" << std::endl;
//		}
}





} // namespace openMVG
} // namespace matching_image_collection
