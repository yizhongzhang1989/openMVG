// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

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
	float distRatio, std::string bin_dir,double maxDistanceThreshold,const sfm::SfM_Data& sfm_data
):Matcher(), f_dist_ratio_(distRatio), bin_dir_(bin_dir),MaxDistanceThreshold(maxDistanceThreshold)
, opticalflow_container(bin_dir, maxDistanceThreshold, sfm_data)
{
	
}

namespace impl
{

// compute matches in local and for every feature in first view ,find least NN features paired in second view.
//      view1											view2
//    -------------------						-------------------
//	  | feat_id1(view1)-+-----------------------+-feat_id2(view2) |
//    |        |        |    feature match		|       |---------+---> distance tracked by optical flow
//    |        ---------+-----------------------+>feat_id1(view1) |
//    -------------------	optical flow		-------------------
//                          
//                         
//
// @param[in]  regions_provider           the information of features.
// @param[in]  pairs                      the view pairs to matching
// @param[in]  fDistRatio                 Distance ratio to discard non meaningful matches
// @param[out] map_PutativesMatches       the computed matches
// @param[in]  opticalflow_container      the optical information 
// @param[in]  my_progress_bar            progress bar ui

// @return                whether the 3d point lies in front of both cameras.
// (Owned by BC)

template <typename ScalarT>
bool Match
(
  const sfm::Regions_Provider & regions_provider,
  const Pair_Set & pairs,
  float fDistRatio,
  PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
  C_Progress * my_progress_bar,
	const OpticalFlow_Container& opticalflow_container
)
{
	
	////statistic////
	std::vector<size_t> vec_local_num_matches;
	std::vector<size_t> vec_num_matches_desfiltering;
	std::vector<size_t> vec_last_num_matches;


  if (!my_progress_bar)
    my_progress_bar = &C_Progress::dummy();
  my_progress_bar->restart(pairs.size(), "\n- Optical Matching -\n");

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

  ///////////////////////////////////////
  ////// Init the cascade hasher  ///////
  ///////////////////////////////////////
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
  ///////////////////////////////////////////////
  //////////////optical matching/////////////////
  ///////////////////////////////////////////////
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

	  //////////////a. find features pairs on whole image(global)//////////////
      IndMatches pvec_indices_global;
      using ResultType = typename Accumulator<ScalarT>::Type;
      std::vector<ResultType> pvec_distances_global;
	  pvec_distances_global.reserve(regionsJ->RegionCount() * 2);
	  pvec_indices_global.reserve(regionsJ->RegionCount() * 2);
	  const std::vector<features::PointFeature> pointFeaturesJ = regionsJ->GetRegionsPositions();
	  // Match the query descriptors to the database
	  cascade_hasher.Match_HashedDescriptions<BaseMat, ResultType>(
		  hashed_base_[J], mat_J,
		  hashed_base_[I], mat_I,
		  &pvec_indices_global, &pvec_distances_global,1);


	  
	  //match by optical flow
	  //std::cout << "Pair " << I << "-" << J << "\n";

	  ////////////b. optical matching /////////////////////
	  //        * find feature pairs on local region limited by tracked distance
	  IndMatches pvec_indices_local;
	  std::vector<ResultType> pvec_distances_local;
	  pvec_distances_local.reserve(regionsJ->RegionCount() * 2);
	  pvec_indices_local.reserve(regionsJ->RegionCount() * 2);
	  bool bofmatching = opticalflow_container.Optical_Matching<BaseMat, ResultType>
		  (
			  J,
			  hashed_base_[J], mat_J,
			  I,
			  hashed_base_[I], mat_I,

			  pointFeaturesJ,
			  pointFeaturesI,
			  (I >= J ? I : J),

			  cascade_hasher.nb_hash_code_,               //cascade_hasher is initialized with the dimension
			  &pvec_indices_local,
			  &pvec_distances_local,
			  1
		 );
	  local_num_matches = pvec_indices_local.size();
	  //
	  ////c. select the feature pairs satisfying ratio threshold between local distance and global distance ////
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

	  std::cout << "The statistic of desfiltering_num_matches:\n";
	  sort(vec_num_matches_desfiltering.begin(), vec_num_matches_desfiltering.end());
	  utils::MinMaxMedianMean<size_t>(std::cout, vec_num_matches_desfiltering);  

	  std::cout << "The statistic of last_num_matches:\n";
	  sort(vec_last_num_matches.begin(), vec_last_num_matches.end());
	  utils::MinMaxMedianMean<size_t>(std::cout, vec_last_num_matches);
	  
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
      f_dist_ratio_,
      map_PutativesMatches,
      my_progress_bar,
	  opticalflow_container);
  }
  else if (regions_provider->Type_id() == typeid(float).name())
  {
	  impl::Match<float>(
		  *regions_provider.get(),
		  pairs,
		  f_dist_ratio_,
		  map_PutativesMatches,
		  my_progress_bar,
		  opticalflow_container);
  }
  else
  {
    std::cerr << "Matcher not implemented for this region type" << std::endl;
  }
}



} // namespace openMVG
} // namespace matching_image_collection
