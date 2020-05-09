// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/sfm/pipelines/global/sfm_global_engine_iterative_pose_augmentation.hpp"
#include "openMVG_IMU/sfm/pipelines/global/myoutput.hpp"
#include "openMVG_IMU/matching_image_collection/ComputeMatchesController.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"

#include "openMVG/sfm/pipelines/relative_pose_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/types.hpp"


#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <ceres/types.h>

#include <iostream>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

GlobalSfMReconstructionEngine_IterativePoseAugmentation::GlobalSfMReconstructionEngine_IterativePoseAugmentation(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : GlobalSfMReconstructionEngine_RelativeMotions_General(sfm_data,soutDirectory,sloggingFile)
{

  
}

GlobalSfMReconstructionEngine_IterativePoseAugmentation::~GlobalSfMReconstructionEngine_IterativePoseAugmentation()
{
  
}

void GlobalSfMReconstructionEngine_IterativePoseAugmentation::SetMatchesDir(const std::string MatchesDir)
{
    sMatchesDir_=MatchesDir;
}

bool GlobalSfMReconstructionEngine_IterativePoseAugmentation::Run()
{
  Pair_Set extra_pairs;
  unsigned int loop_i = 1;
  //Because the matching of a image pair is only based on the features in the image,
  //even though optimized by ba, the failed image pair in matching phase will still fail.
  //So we develop a set to record the failed or succeeded matches to avoid duplicate matching and boost speeds
  Pair_Set tried_pairs;
  while(LoopDetection(extra_pairs, tried_pairs))
  {
	  system::Timer augmentation_once_timer;
      std::cout<<"///////////////////////////\n"
               <<"//////////loop "<<loop_i<<"/////////\n"
               <<"///////////////////////////\n";
	  tried_pairs.insert(extra_pairs.begin(), extra_pairs.end());
      std::cout<<"Detect "<<extra_pairs.size()<<" pairs\n";
      std::shared_ptr<Matches_Provider> extra_matches_provider = std::make_shared<Matches_Provider>();
	 
	  system::Timer matching_timer;
	  matching_image_collection::ComputeMatchesController::Process(sfm_data_,sMatchesDir_,extra_matches_provider.get(),"f",
                                       true,extra_pairs,true, stlplus::create_filespec(sOut_directory_, "matches_f_" + std::to_string(loop_i), ".bin"));
	  std::cout << "Matching task done in (s): " << matching_timer.elapsed() << std::endl;

      std::cout<<"Compute "<<extra_matches_provider->pairWise_matches_.size()<<" matches in "<<extra_pairs.size()<<" pairs\n";
      if(extra_matches_provider->pairWise_matches_.size()==0)
      {
        std::cout<<"No more matches are found\n";
        break;
      }
      // add new extra matches into total matches
      matching::PairWiseMatches::iterator iter;
      std::cout << "matches_provider->getpairs() first" << matches_provider_->pairWise_matches_.size() << "\n";
      std::cout << "extra_matches_provider->getpairs() first" << extra_matches_provider->pairWise_matches_.size() << "\n";
      for (iter = extra_matches_provider->pairWise_matches_.begin(); iter != extra_matches_provider->pairWise_matches_.end();
        )
      {
        if (matches_provider_->pairWise_matches_.count(iter->first) == 0)
        {
          matches_provider_->pairWise_matches_.insert(*iter);
          iter++;
        }
        else  //remove already existing matches 
        {
          iter=extra_matches_provider->pairWise_matches_.erase(iter);
        }
      }
      //matches_provider now contains all matches(original and extra)
      std::cout << "matches_provider->getpairs() then" << matches_provider_->pairWise_matches_.size() << "\n";
      std::cout << "extra_matches_provider->getpairs() then" << extra_matches_provider->pairWise_matches_.size() << "\n";

      if(extra_matches_provider->pairWise_matches_.size()==0)
      {
        std::cout<<"No more matches are found\n";
        break;
      }
	  ////BC start////
	  //save the raw pose
	  Output_trajectory(stlplus::create_filespec(sOut_directory_, "rawposes_" + std::to_string(loop_i), ".csv"), sfm_data_);

	  Output_Matchings(stlplus::create_filespec(sOut_directory_, "matching_withloopbegin_" + std::to_string(loop_i), ".csv"), matches_provider_);
	  Output_Matchings(stlplus::create_filespec(sOut_directory_, "extramatching_withloopbegin_" + std::to_string(loop_i), ".csv"), extra_matches_provider.get());
	  //save the raw corresponds

	  ////BC end //// 
	  system::Timer averaging_timer;
      if(!Process(extra_matches_provider))
      {
        std::cout<<"Global sfm augmentation failed\n";
        break;
      }
	  std::cout << "Remain " << extra_matches_provider->pairWise_matches_.size() << " extra pairs\n";
	  std::cout << "Averaging task done in (s): " << averaging_timer.elapsed() << std::endl;
	  /////////////BC start////////////////

	  //save the image pose by view_id
	  Output_trajectory(stlplus::create_filespec(sOut_directory_, "view_poses_" + std::to_string(loop_i), ".csv"), sfm_data_);

	  //save the triangulated correspondings
	  Output_TriangulatedCorrespondings(stlplus::create_filespec(sOut_directory_, "triangulated_correspondings_" + std::to_string(loop_i), ".csv"), sfm_data_);

	  // save triangulated matchings
	  Output_TriangulatedMatchings(stlplus::create_filespec(sOut_directory_, "triangulated_matchings_" + std::to_string(loop_i), ".csv"), matches_provider_, sfm_data_);
	  //save remain extra matches
	  Output_Matchings(stlplus::create_filespec(sOut_directory_, "remain_extramatchings_" + std::to_string(loop_i), ".csv"), extra_matches_provider.get());

	  /////////////BC start////////////////
      //optimize
	  system::Timer optimizing_timer;
	  Optimize();
	  std::cout << "Optimizing task done in (s): " << optimizing_timer.elapsed() << std::endl;
	  Output_trajectory(stlplus::create_filespec(sOut_directory_, "view_poses_ba_" + std::to_string(loop_i), ".csv"), sfm_data_);

	  Output_TriangulatedCorrespondings(stlplus::create_filespec(sOut_directory_, "triangulated_correspondings_ba_" + std::to_string(loop_i), ".csv"), sfm_data_);
      loop_i++;
	  std::cout << "This augmentation task done in (s): " << augmentation_once_timer.elapsed() << std::endl;

  }
  return true; 
}

bool GlobalSfMReconstructionEngine_IterativePoseAugmentation::LoopDetection(Pair_Set& extra_pairs,const Pair_Set& tried_pairs)
{
  extra_pairs.clear();
	const double MaxAngleThreshold = 25.0;
  for(const auto& view_i:sfm_data_.GetViews())
  {
    const Pose3& pose_i = sfm_data_.GetPoseOrDie(view_i.second.get());
    const Mat3 R_i = pose_i.rotation();
    for(const auto& view_j:sfm_data_.GetViews())
    {
      if(view_i.first>=view_j.first) continue;
      Pair pair(view_i.first,view_j.first);
	  Pair pair_inverse(view_j.first, view_i.first);
      if(matches_provider_->pairWise_matches_.count(pair)||
		 matches_provider_->pairWise_matches_.count(pair_inverse)) continue;
	  if (tried_pairs.count(pair) || tried_pairs.count(pair_inverse)) continue;

      const Pose3& pose_j = sfm_data_.GetPoseOrDie(view_j.second.get());
      const Mat3 R_j = pose_j.rotation();
      const double angularErrorDegree = R2D(getRotationMagnitude(R_i * R_j.transpose()));
      if (angularErrorDegree < MaxAngleThreshold)
      {
          extra_pairs.insert(pair);
      }
    }
  }
  return extra_pairs.size() > 0;
}

bool GlobalSfMReconstructionEngine_IterativePoseAugmentation::Optimize() {

   std::cout<<"/////IMU Global SfM BA/////\n";
  
  if (!Adjust())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }
  
  //-- Export statistics about the SfM process
  std::cout << "Structure from Motion statistics.";


  
  std::cout << "-------------------------------" << "\n"
    << "-- View count: " << sfm_data_.GetViews().size() << "\n"
    << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "\n"
    << "-- Pose count: " << sfm_data_.GetPoses().size() << "\n"
    << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";
    

  return true;
}

bool GlobalSfMReconstructionEngine_IterativePoseAugmentation::Process(std::shared_ptr<Matches_Provider> extra_matches_provider_) {


  
  std::cout<<"/////IMU Global SfM Iterative Pose Augmentation/////\n";
  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }
  ////bc start////
  if(sfm_data_.GetPoses().size()==0)
  {
    std::cout<<"Invalid input pose graph for global sfm\n";
    return false;
  }

  Hash_Map<IndexT, Mat3> global_rotations;
  for(const auto& pose_item: sfm_data_.GetPoses())
  {
      global_rotations.emplace(pose_item.first, pose_item.second.rotation());
  }
  ////bc end////
  matching::PairWiseMatches  tripletWise_matches;
  if (!Compute_Global_PrioTranslations(global_rotations, tripletWise_matches,extra_matches_provider_.get()))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  //recompute structure
  sfm_data_.structure.clear();    //bc
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  

  return true;
}


/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_IterativePoseAugmentation::Compute_Global_PrioTranslations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches,
  Matches_Provider* extra_matches_provider_
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_PrioTranslation_AveragingSolver priotranslation_averaging_solver;
  const bool bTranslationAveraging = priotranslation_averaging_solver.MyRun(
    eTranslation_averaging_method_,
    sfm_data_,
    features_provider_,
    matches_provider_,
    global_rotations,
    tripletWise_matches,
    extra_matches_provider_);
  /////////////BC start////////////////
  //output the pose graph
  std::cout<<"BC is debuging\n";
  std::ofstream extratripletestimation_file(stlplus::create_filespec(sOut_directory_, "extratriplet_estimation.csv"));
  for (const auto& pair_item : priotranslation_averaging_solver.extra_pairs_)
  {
	  extratripletestimation_file << pair_item.first.first << "-" << pair_item.first.second << ":" << pair_item.second << "\n";
  }
  extratripletestimation_file.close();
  if (!sLogging_file_.empty() && !sOut_directory_.empty())
  {
    std::cout<<"saving pose graph after translation\n";
    std::ofstream viewgraph_file(stlplus::create_filespec(sOut_directory_, "posegraph_augmentation.csv"));
    std::set<IndexT> set_view_ids;
    Pair_Set relative_view_pairs;

    viewgraph_file<<"pose_id1,pose_id2\n";
    
      for (const openMVG::RelativeInfo_Vec & iter : priotranslation_averaging_solver.Getrelative_motion())
      {
          for (const relativeInfo & rel : iter)
          {
            relative_view_pairs.insert(rel.first);
            set_view_ids.insert(rel.first.first);
            set_view_ids.insert(rel.first.second);
            viewgraph_file<<rel.first.first<<","<<rel.first.second<<"\n";
          }
      }
    
    // Log a relative view graph
    {
      const std::string sGraph_name = "posegraph_augmentation";
      graph::indexedGraph putativeGraph(set_view_ids, relative_view_pairs);
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, sGraph_name),
        putativeGraph);
    }
    viewgraph_file.close();
  }
  std::cout<<"Debug complete\n";
  /////////////BC   end////////////////

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}



} // namespace sfm
} // namespace openMVG
