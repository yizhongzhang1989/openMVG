// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG_IMU/sfm/pipelines/global/sfm_global_engine_pose_augmentation.hpp"

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

GlobalSfMReconstructionEngine_PoseAugmentation::GlobalSfMReconstructionEngine_PoseAugmentation(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : GlobalSfMReconstructionEngine_RelativeMotions_General(sfm_data,soutDirectory,sloggingFile)
{

  
}

GlobalSfMReconstructionEngine_PoseAugmentation::~GlobalSfMReconstructionEngine_PoseAugmentation()
{
  
}


////bc start/////
void GlobalSfMReconstructionEngine_PoseAugmentation::SetExtraMatchesProvider(Matches_Provider * provider)
{
  extra_matches_provider_ = provider;
}
////bc end/////


//int GlobalSfMReconstructionEngine_PoseAugmentation::RemainExtraPars(const auto& )

bool GlobalSfMReconstructionEngine_PoseAugmentation::Process() {

   std::cout<<"/////IMU Global SfM Pose Augmentation/////\n";
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
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  /////////////BC start////////////////
  //output the view  graph
  std::cout<<"BC is debuging\n";
  if (!sLogging_file_.empty() && !sOut_directory_.empty())
  {
    std::cout<<"saving view graph before final ba\n";
    std::ofstream viewgraph_file(stlplus::create_filespec(sOut_directory_, "viewgraph_augamentation.csv"));
    std::set<IndexT> set_view_ids;
    Pair_Set relative_view_pairs;

    viewgraph_file<<"image_id1,image_id2,match_num\n";
    for(const auto& pair_iter : tripletWise_matches)
    {
      const auto & imageI = pair_iter.first.first;
      const auto & imageJ = pair_iter.first.second;
      size_t match_num = pair_iter.second.size();
      set_view_ids.insert(imageI);
      set_view_ids.insert(imageJ);
      relative_view_pairs.insert(Pair(imageI,imageJ));

      viewgraph_file<<imageI<<","<<imageJ<<","<<match_num<<"\n";
    } 
    // Log a relative view graph
    {
      const std::string sGraph_name = "viewgraph_augmentation";
      graph::indexedGraph putativeGraph(set_view_ids, relative_view_pairs);
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, sGraph_name),
        putativeGraph);
    }
	viewgraph_file.close();
  }
  //save the sfm scene data
  Save(sfm_data_,
	  stlplus::create_filespec(sOut_directory_, "Trajectory_Augmentation", ".json"),
	  ESfM_Data(INTRINSICS | VIEWS | EXTRINSICS));
  Save(sfm_data_,
      stlplus::create_filespec(sOut_directory_, "sfm_data_Augmentation",".json"),
      ESfM_Data(ALL));
  
  std::cout<<"Debug complete\n";
  /////////////BC   end////////////////
  // if (!Adjust())
  // {
  //   std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
  //   return false;
  // }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}


/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_PoseAugmentation::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches
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
