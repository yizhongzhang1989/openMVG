// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

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
	bool bfeature_validation, bool bopticalfiltering,
    bool bdynamicdistance, bool bopticalmatching,bool bnotablefeaturesdetection,bool bnotablevalidation
):Matcher(), f_dist_ratio_(distRatio),
opticalflow_container(bin_dir,maxDistanceThreshold,sfm_data),
bfeature_validation_(bfeature_validation), bopticalfiltering_(bopticalfiltering),
bdynamicdistance_(bdynamicdistance), bopticalmatching_(bopticalmatching),
bnotablefeaturesdetection_(bnotablefeaturesdetection), bnotablevalidation_(bnotablevalidation)
{
	
	
}

namespace impl
{

// hierarchical matching
// feature matching// *	sift matching// *	validation : I->J && J->I    <- bfeature_validation// optical filtering                 <- bopticalfiltering// *	static threshold // *    dynamic threshold            <- bdynamicdistance
// notable features                  <- bnotablefeaturesdetection
// *	compute notable features
// *	validation                   <- bnotablevalidation
// optical matching                  <- bopticalmatching
// *	local matching
// *	validation
//                          
//
// @param[in]  regions_provider           the information of features.
// @param[in]  pairs                      the view pairs to matching
// @param[in]  fDistRatio                 Distance ratio to discard non meaningful matches
// @param[out] map_PutativesMatches       the computed matches
// @param[in]  opticalflow_container      the optical information 
// @param[in]  my_progress_bar            progress bar ui
// @param[out] notable_features           computed notable features
// @param[in]  bfeature_validation        signal of the duplex validation of matches from I to J and J to I
// @param[in]  bopticalfiltering          signal of optical filtering
// @param[in]  bdynamicdistance           signal of using dynamic threshold
// @param[in]  bopticalmatching           signal of optical matching
// @param[in]  bnotablefeaturesdetection  signal of computing notable features
// @param[in]  bnotablevalidation         signal of notable validation
// (Owned by BC)

template <typename ScalarT>
void Match_Hierarchical
(
	const sfm::Regions_Provider & regions_provider,
	const Pair_Set & pairs,
	float fDistRatio,
	PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
	OpticalFlow_Container& opticalflow_container,
	C_Progress * my_progress_bar,
	std::map<IndexT,std::set<IndexT>>& notable_features,
	bool bfeature_validation = true,
	bool bopticalfiltering = true,
	bool bdynamicdistance = false,
	bool bopticalmatching = false,
	bool bnotablefeaturesdetection = false,
	bool bnotablevalidation = false
)
{
  
  using ResultType = typename Accumulator<ScalarT>::Type;
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  ////statistic////
  std::vector<size_t> first_num_matches;
  std::vector<size_t> second_num_matches;
  std::vector<size_t> third_num_matches;
  std::vector<size_t> fourth_num_matches;
  std::vector<size_t> fifth_num_matches;
  std::vector<size_t> sixth_num_matches;
  std::vector<size_t> seventh_num_matches;
  std::vector<size_t> eighth_num_matches;
  //std::map<Pair, std::map<Pair, ResultType>> descdis_putativematches;
  //END(Author: BC)===================================================

  
  if (!my_progress_bar)
	  my_progress_bar = &C_Progress::dummy();
  my_progress_bar->restart(pairs.size(), "\n- Matching -\n");
  
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
  ///////////////////////////////////////
  ////// Init the cascade hasher  ///////
  ///////////////////////////////////////
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  CascadeHasher2 cascade_hasher2;
  //END(Author: BC)===================================================
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
  ///////////////////////////////////////
  ////////// Sift matching  /////////////
  ///////////////////////////////////////
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
	
	
	//////// feature validation ///////////
	////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
	//find one feature in J for every feature in I
	if(bfeature_validation)
	{
		std::set<IndexT> matched_featIids;
		std::vector<ResultType> pvec_distances_IJ;
		// get the matches in first matching(J->I)
		for (size_t k=0; k < vec_nn_ratio_idx_JI.size(); ++k)
		{
			const size_t index = vec_nn_ratio_idx_JI[k];
			matched_featIids.insert(pvec_indices_JI[index*2].j_);
		}
		pvec_distances_IJ.reserve(matched_featIids.size() * 2);
		pvec_indices_IJ.reserve(matched_featIids.size() * 2);
		//only find matches for first matches pair from I to J
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

	//check the feautre pair both appear in two kinds of matches.
	matching::IndMatches vec_putative_matches;
	std::map<Pair, ResultType> map_putative_descdis;
	bool bindex_consistent = true;  //check the correctness of features
	vec_putative_matches.reserve(vec_nn_ratio_idx_JI.size());
	for (size_t k = 0; k < vec_nn_ratio_idx_JI.size(); ++k)
	{
		const size_t index1 = vec_nn_ratio_idx_JI[k];
		const IndexT feat_i_id1 = pvec_indices_JI[index1 * 2].j_;
		const IndexT feat_j_id1 = pvec_indices_JI[index1 * 2].i_;
		const features::PointFeature pointfeati_1 = pointFeaturesI[feat_i_id1];
		const features::PointFeature pointfeatj_1 = pointFeaturesJ[feat_j_id1];
		bool boccur_featurepair = false;
		if (!bfeature_validation)  // accept the feature pair if user does not specify feature validation
		{
			boccur_featurepair = true;
		}
		else
		{
			//check whether the feature pair is found from I to J
			for (size_t w = 0; w < vec_nn_ratio_idx_IJ.size(); ++w)
			{
				const size_t index2 = vec_nn_ratio_idx_IJ[w];
				const IndexT feat_i_id2 = pvec_indices_IJ[index2 * 2].i_;
				const features::PointFeature pointfeati_2 = pointFeaturesI[feat_i_id2];
				if (feat_i_id2 != feat_i_id1 && pointfeati_2.coords() != pointfeati_1.coords()) continue;

				const IndexT feat_j_id2 = pvec_indices_IJ[index2 * 2].j_;
				
				const features::PointFeature pointfeatj_2 = pointFeaturesJ[feat_j_id2];
				//check whether the feature in view J is consistent
				if (feat_j_id1 == feat_j_id2 || (pointfeatj_1.coords() == pointfeatj_2.coords()))
				{
					boccur_featurepair = true;
					break;
				}
				//break;
			}
		}
		//the pair are found both from I to J and from J to I
		if (boccur_featurepair)
		{
			vec_putative_matches.emplace_back(feat_i_id1, feat_j_id1);
			map_putative_descdis.emplace(Pair(feat_i_id1, feat_j_id1), pvec_distances_JI[index1 * 2]);
		}

	}
	num_matches_second = vec_putative_matches.size();
	//END(Author: BC)===================================================

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
		if (!vec_putative_matches.empty())
		{

			map_PutativesMatches.insert(
			{
			  {I,J},
			  std::move(vec_putative_matches)
			});
			////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
			//descdis_putativematches.emplace(Pair(I, J), std::move(map_putative_descdis));
			first_num_matches.push_back(num_matches_first);
			second_num_matches.push_back(num_matches_second);
			//END(Author: BC)===================================================
		}
	}
	++(*my_progress_bar);
	}
  }
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  //Display some logs
  std::cout << "The statistic of first_num_matches:\n";
  sort(first_num_matches.begin(), first_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, first_num_matches);
  std::cout << "The statistic of second_num_matches:\n";
  sort(second_num_matches.begin(), second_num_matches.end());
  utils::MinMaxMedianMean<size_t>(std::cout, second_num_matches);
  
	////////////////////////////////////////////////
	//////////////optical filtering/////////////////
	////////////////////////////////////////////////
	//check whether the optical data are read successfully
	if (!opticalflow_container.readable_)
	{
		  std::cout << "the binary files are not read completely.\n";
		  return;
	}
	// if bdynamicdistance is true , set dynamic threshold
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
	if (bopticalfiltering)
	{

			// optical filtering
			bool bfiltering = opticalflow_container.Optical_Filtering
			(map_PutativesMatches, &regions_provider, my_progress_bar
			);

			if (!bfiltering)
			{
				return;
			}
			/// if bnotablefeaturesdetection is true, compute notable features
			if(bnotablefeaturesdetection)
			{
				/////compute notable features in every frame that has matches with next 2 frame after optical filtering.
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
					
				}

			}
			//{
			//	//delete the distance records of filtered feature pair and update ``
			//	for (const auto& pairwise_item : map_PutativesMatches)
			//	{

			//		std::set<Pair> existed_featpair;
			//		for (const auto& indmatch : pairwise_item.second)
			//		{
			//			existed_featpair.insert(Pair(indmatch.i_, indmatch.j_));
			//		}
			//		std::map<Pair, ResultType>& map_descdis_IJ = descdis_putativematches.at(pairwise_item.first);
			//		for (std::map<Pair, ResultType>::iterator it = map_descdis_IJ.begin(); it != map_descdis_IJ.end();)
			//		{
			//			if (!existed_featpair.count(it->first))
			//			{
			//				it = map_descdis_IJ.erase(it);
			//			}
			//			else
			//			{
			//				it++;
			//			}
			//		}

			//	}
			//}
			
			{
				std::cout << "The statistic of third_num_matches:\n";
				sort(third_num_matches.begin(), third_num_matches.end());
				utils::MinMaxMedianMean<size_t>(std::cout, third_num_matches);
				std::cout << "The statistic of fourth_num_matches:\n";
				sort(fourth_num_matches.begin(), fourth_num_matches.end());
				utils::MinMaxMedianMean<size_t>(std::cout, fourth_num_matches);
			}
	}
	  ///////////////////////////////////////////////
	  //////////////optical matching/////////////////
	  ///////////////////////////////////////////////
	if (bopticalmatching)
	{


		std::cout << "optical threshold: " << opticalflow_container.MaxDistanceThreshold_ << "\n";
		typedef std::pair<IndexT, std::vector<IndMatch>*> Pairiva;
		std::map<IndexT, std::vector<Pairiva>> map_putmatches;
		//adjust the format of matches for OPENMP
		for (auto& putmatch_item : map_PutativesMatches)
		{
			map_putmatches[putmatch_item.first.first].emplace_back(
				Pairiva(putmatch_item.first.second, &putmatch_item.second));
		}
		if (!my_progress_bar)
			my_progress_bar = &C_Progress::dummy();
		my_progress_bar->restart(map_PutativesMatches.size(), "\n- Optical Matching -\n");
		
		for (auto& putmatch_item : map_putmatches)
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
			for (int j = 0; j < (int)putmatch_item.second.size(); j++)
			{
				if (my_progress_bar->hasBeenCanceled())
					continue;



				Pairiva& pairiva = putmatch_item.second[j];
				const IndexT J = pairiva.first;
				//////// find dynamic descriptor similarity threshold for each image pair /////
				////       * select the maximum descriptor distance of existing matches as thershold
				//ResultType MaxSimilarity_perpair = 0.0;
				//for (const auto& descdis_item : descdis_putativematches.at(Pair(I, J)))
				//{
				//	MaxSimilarity_perpair = MaxSimilarity_perpair > descdis_item.second ?
				//		MaxSimilarity_perpair : descdis_item.second;
				//}

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
				//////////////a. find features pairs on whole image(global)//////////////
				IndMatches pvec_indices_global;

				std::vector<ResultType> pvec_distances_global;
				pvec_distances_global.reserve(regionsJ->RegionCount() * 2);
				pvec_indices_global.reserve(regionsJ->RegionCount() * 2);

				cascade_hasher2.Match_HashedDescriptions<BaseMat, ResultType>(
					hashed_base_[J], mat_J,
					hashed_base_[I], mat_I,
					&pvec_indices_global, &pvec_distances_global, 1);
				//record the features that have been paired.
				for (const auto& featmatch_item : (*pairiva.second))
				{
					paired_featiid.insert(featmatch_item.i_);
					paired_featjid.insert(featmatch_item.j_);
				}
				std::vector<ResultType> pvec_distances_JI;
				IndMatches pvec_indices_JI;
				////////////b. optical matching /////////////////////
				//        * find feature pairs on local region limited by tracked distance

				bool bofmatching = opticalflow_container.Optical_Matching<BaseMat, ResultType>
					(
						J,
						hashed_base_[J], mat_J,
						I,
						hashed_base_[I], mat_I,

						pointFeaturesJ,
						pointFeaturesI,
						
						(I >= J ? I : J),

						cascade_hasher2.nb_hash_code_,               //cascade_hasher is initialized with the dimension
						&pvec_indices_JI,
						&pvec_distances_JI,
						
						1,
						&paired_featjid,
						&paired_featiid
						);
				num_fifth_matches = pvec_indices_JI.size() + num_fourth_matches;
				matching::IndMatches vec_putative_matches;
				// set threshold  as MaxSimilarity_perpair
				/*for (size_t k = 0; k < pvec_indices_JI.size(); ++k)
				{
					if (pvec_distances_JI[k] < MaxSimilarity_perpair)
						vec_putative_matches.emplace_back(pvec_indices_JI[k].j_, pvec_indices_JI[k].i_);
				}*/

				////c. select the feature pairs satisfying ratio threshold between local distance and global distance ////
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

					pairiva.second->insert(pairiva.second->end(), vec_putative_matches.begin(), vec_putative_matches.end());
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
  //END(Author: BC)===================================================
  
}
} // namespace impl



void Hierarchical_Matcher_Regions::Match_Hierarchical
(
	const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
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

	if (regions_provider->Type_id() == typeid(unsigned char).name())
	{
		impl::Match_Hierarchical<unsigned char>(
			*regions_provider.get(),
			pairs,
			f_dist_ratio_,
			map_PutativesMatches,
			opticalflow_container,
			my_progress_bar, 
			notable_features,
			bfeature_validation_,
			bopticalfiltering_,
			bdynamicdistance_,
			bopticalmatching_, bnotablefeaturesdetection_, bnotablevalidation_);
	}
	else if (regions_provider->Type_id() == typeid(float).name())
	{
		impl::Match_Hierarchical<float>(
			*regions_provider.get(),
			pairs,
			f_dist_ratio_,
			map_PutativesMatches,
			opticalflow_container,
			my_progress_bar,
			notable_features, 
			bfeature_validation_,
			bopticalfiltering_,
			bdynamicdistance_,
			bopticalmatching_, bnotablefeaturesdetection_, bnotablevalidation_);
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

}





} // namespace openMVG
} // namespace matching_image_collection
