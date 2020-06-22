// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_stats.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust_Angular.hpp"
#include "openMVG/matching_image_collection/Eo_Robust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "openMVG_IMU/matching_image_collection/Optical_Flow_Matcher_Regions.hpp"
#include "openMVG_IMU/matching_image_collection/Hierarchical_Matcher_Regions.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <queue>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;
using namespace std;

enum EGeometricModel
{
  FUNDAMENTAL_MATRIX = 0,
  ESSENTIAL_MATRIX   = 1,
  HOMOGRAPHY_MATRIX  = 2,
  ESSENTIAL_MATRIX_ANGULAR = 3,
  ESSENTIAL_MATRIX_ORTHO = 4
};

enum EPairMode
{
  PAIR_EXHAUSTIVE = 0,
  PAIR_CONTIGUOUS = 1,
  PAIR_FROM_FILE  = 2
};

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
/// - Export computed data

void SolveArticulationPoint(const PairWiseMatches& map_GeometricMatches, const int i_neighbour_count
	, Pair_Set& pairs)
{


	//collect vertex
	std::set<IndexT> used_index;
	for (const auto& pairwisematches : map_GeometricMatches)
	{
		used_index.insert(pairwisematches.first.first);
		used_index.insert(pairwisematches.first.second);
	}

	//Discretization
	std::map<IndexT, IndexT> reindex_front;
	std::map<IndexT, IndexT> reindex_back;
	size_t num_view = 0;
	for (const IndexT view_id : used_index)
	{
		reindex_front.emplace(view_id, num_view);
		reindex_back.emplace(num_view, view_id);
		num_view++;
	}

	//collect edges
	std::vector<std::vector<IndexT>> edges;
	edges.resize(num_view);
	for (const auto& pairwisematches : map_GeometricMatches)
	{
		IndexT u = reindex_front.at(pairwisematches.first.first);
		IndexT v = reindex_front.at(pairwisematches.first.second);
		edges[u].push_back(v);
		edges[v].push_back(u);
	}
	//travel all node by depth first search
	auto dfs_function = [&](IndexT i, std::vector<bool>& vis, const std::vector<std::vector<IndexT>>& edges)
	{
		std::stack<IndexT> dfs_stack;
		dfs_stack.push(i);
		vis[i] = true;
		while (!dfs_stack.empty())
		{
			IndexT u = dfs_stack.top();
			dfs_stack.pop();
			for (const auto& v : edges[u])
			{
				if (vis[v]) continue;
				dfs_stack.push(v);
				vis[v] = true;
			}
		}
	};
	//find number of connected components in graph
	size_t num_connected = 0;
	std::vector<bool> vis(num_view);
	std::fill(vis.begin(), vis.end(), false);
	for (size_t i = 0; i < num_view; i++)
	{
		if (!vis[i])
		{
			dfs_function(i, vis, edges);
			num_connected++;
		}
	}
	std::cout << "#biconnected componcents: " << num_connected << "\n";
	//find articulation point in graph
	std::map<IndexT, int> articulation_points;

	for (size_t i = 0; i < num_view; i++)
	{
		std::fill(vis.begin(), vis.end(), false);
		size_t num_conected_tmp = 0;
		vis[i] = true;  //remove the point 
						//check whether new connected components occur 
		for (size_t j = 0; j < num_view; j++)
		{
			if (!vis[j])
			{
				dfs_function(j, vis, edges);
				num_conected_tmp++;
			}
		}
		if (num_conected_tmp != num_connected)
		{
			articulation_points.emplace(i, num_conected_tmp -num_connected );
			std::cout << "find articulation points: view " << reindex_back.at(i) << "\n";
		}
	}


	//   *** -- * -- ***  middle point is articulation point
	// we should add more edges between left nodes and right nodes.
	//xxxxxPlan1 : pair the articulation views with its neighbour frames (optional)
	//Plan2 : pair the articulation views according graph
	int num_neighbour = i_neighbour_count < 0 ? 15 : i_neighbour_count;



	for (const auto& node_item : articulation_points)
	{

		IndexT node_id = node_item.first;
		IndexT view_id = reindex_back.at(node_id);
		//IndexT first_view = view_id > i_neighbour_count ? view_id - i_neighbour_count : 0;
		//IndexT second_view = view_id + i_neighbour_count < *used_index.rbegin() ? view_id + i_neighbour_count : *used_index.rbegin();
		std::set<IndexT> searched_views;
		searched_views.insert(view_id);
		int minSearchDepth = num_neighbour + 1;
		int MinSearchNodes = num_neighbour * (node_item.second + 1);
		//std::cout << "MinSearchNodes:" << MinSearchNodes << "\n";
		{
			std::queue<IndexT> bfs_queue;
			std::vector<int> vis_tmp(num_view);
			std::fill(vis_tmp.begin(), vis_tmp.end(), 0);
			bfs_queue.push(node_id);
			vis_tmp[node_id] = 1;
			while (!bfs_queue.empty())
			{
				IndexT u = bfs_queue.front();
				bfs_queue.pop();
				for (const auto& v : edges[u])
				{
					if (vis_tmp[v]) continue;

					vis_tmp[v] = vis_tmp[u] + 1;
					if (vis_tmp[v] <= minSearchDepth)
					{
						bfs_queue.push(v);
						//std::cout << "find view" << reindex_back.at(v) << "in level " << vis_tmp[v] << "\n";
						//std::cout << "			minSearchDepth:" << minSearchDepth << "\n";
						searched_views.insert(reindex_back.at(v));
						if (searched_views.size() > MinSearchNodes)
						{
							minSearchDepth = vis_tmp[v];
						}
					}

				}
			}
		}
		size_t  size_before = pairs.size();
		for (const auto& view_id1 : searched_views)
		{
			for (const auto& view_id2 : searched_views)
			{
				if (view_id1 < view_id2)
				{
					pairs.insert(Pair(view_id1, view_id2));
				}
			}
		}
		size_t  size_after = pairs.size();
		std::cout << "find " << size_after - size_before << " pairs for articulation view " << view_id << "\n";

	}
	std::cout << "#total pairs found : " << pairs.size() << "\n";
}



int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDirectory = "";
  std::string sGeometricModel = "f";
  float fDistRatio = 0.8f;
  int iMatchingVideoMode = -1;
  std::string sPredefinedPairList = "";
  std::string sNearestMatchingMethod = "AUTO";
  std::string bin_dir = "";   //bc
  bool bForce = false;
  bool bGuided_matching = false;
  int imax_iteration = 2048;
  unsigned int ui_max_cache_size = 0;
  double MaxDistanceThreshold = 10.0;
  bool bfeature_validation = true, bopticalfiltering = true;
  bool bdynamicdistance = false, bopticalmatching = false;
  bool bSolveArticulationPoint = true;

  //required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sMatchesDirectory, "out_dir") );
  // Options
  cmd.add( make_option('r', fDistRatio, "ratio") );
  cmd.add( make_option('g', sGeometricModel, "geometric_model") );
  cmd.add( make_option('v', iMatchingVideoMode, "video_mode_matching") );
  cmd.add( make_option('l', sPredefinedPairList, "pair_list") );
  cmd.add( make_option('n', sNearestMatchingMethod, "nearest_matching_method") );
  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('m', bGuided_matching, "guided_matching") );
  cmd.add( make_option('I', imax_iteration, "max_iteration") );
  cmd.add( make_option('c', ui_max_cache_size, "cache_size") );
  cmd.add(make_option('b', bin_dir, "bin_dir"));
  cmd.add(make_option('k', MaxDistanceThreshold, "maxdistancethreshold"));
  cmd.add(make_option('a', bfeature_validation, "bfeature_validation"));
  cmd.add(make_option('d', bopticalfiltering, "bopticalfiltering"));
  cmd.add(make_option('e', bdynamicdistance, "bdynamicdistance"));
  cmd.add(make_option('h', bopticalmatching, "bopticalmatching"));
  cmd.add(make_option('j', bSolveArticulationPoint, "bSolveArticulationPoint"));


  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file\n"
      << "[-o|--out_dir path] output path where computed are stored\n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data]\n"
      << "[-r|--ratio] Distance ratio to discard non meaningful matches\n"
      << "   0.8: (default).\n"
      << "[-g|--geometric_model]\n"
      << "  (pairwise correspondences filtering thanks to robust model estimation):\n"
      << "   f: (default) fundamental matrix,\n"
      << "   e: essential matrix,\n"
      << "   h: homography matrix.\n"
      << "   a: essential matrix with an angular parametrization,\n"
      << "   o: orthographic essential matrix.\n"
      << "[-v|--video_mode_matching]\n"
      << "  (sequence matching with an overlap of X images)\n"
      << "   X: with match 0 with (1->X), ...]\n"
      << "   2: will match 0 with (1,2), 1 with (2,3), ...\n"
      << "   3: will match 0 with (1,2,3), 1 with (2,3,4), ...\n"
      << "[-l]--pair_list] file\n"
      << "[-n|--nearest_matching_method]\n"
      << "  AUTO: auto choice from regions type,\n"
      << "  For Scalar based regions descriptor:\n"
      << "    BRUTEFORCEL2: L2 BruteForce matching,\n"
      << "    HNSWL2: L2 Approximate Matching with Hierarchical Navigable Small World graphs,\n"
      << "    ANNL2: L2 Approximate Nearest Neighbor matching,\n"
      << "    CASCADEHASHINGL2: L2 Cascade Hashing matching.\n"
      << "    FASTCASCADEHASHINGL2: (default)\n"
      << "      L2 Cascade Hashing with precomputed hashed regions\n"
      << "     (faster than CASCADEHASHINGL2 but use more memory).\n"
	  << "     OPTICALFLOW: optical matching\n"
	  << "     HIERARCHICAL: sift matching (+feature validation) (+static/dynamic optical filtering)\n"
	  << "                   (+optical matching).\n"
      << "  For Binary based descriptor:\n"
      << "    BRUTEFORCEHAMMING: BruteForce Hamming matching.\n"
      << "[-m|--guided_matching]\n"
      << "  use the found model to improve the pairwise correspondences.\n"
      << "[-c|--cache_size]\n"
      << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
      << "  If not used, all regions will be load in memory."
      << "[-b|--bin_dir]\n"
	  << "  the directory stores the binary file of position of features tracked in optical flow.\n"
	  << "  It must be specified if you set nearest matching method as `OPTICALFLOW`.\n"
    << "[-k|--maxdistancethreshold]\n"
	    << "  the distance threshold for the optical filtering.\n"
		  << "[-a|--bfeature_validation]\n"
		  << "  matching with I to J and J to I.\n"
		  << "[-d|--bopticalfiltering]\n"
		  << "  filter matches by optical flow.\n"
		  << "[-e|--bdynamicdistance]\n"
		  << "  select a dynamic distance threshold in optical filtering.\n"
		  << "[-h|--bopticalmatching]\n"
		  << "  match with local feature according optical flow.\n"
		  << "[-j|--bSolveArticulationPoint]\n"
		  << "  detect articulation point and repair the graph.\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " << "\n"
			<< argv[0] << "\n"
			<< "--input_file " << sSfM_Data_Filename << "\n"
			<< "--out_dir " << sMatchesDirectory << "\n"
			<< "Optional parameters:" << "\n"
			<< "--force " << bForce << "\n"
			<< "--ratio " << fDistRatio << "\n"
			<< "--geometric_model " << sGeometricModel << "\n"
			<< "--video_mode_matching " << iMatchingVideoMode << "\n"
			<< "--pair_list " << sPredefinedPairList << "\n"
			<< "--nearest_matching_method " << sNearestMatchingMethod << "\n"
			<< "--guided_matching " << bGuided_matching << "\n"
			<< "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << std::endl
			<< "--bin_dir " << bin_dir << "\n"
			<< "--MaxDistanceThreshold " << MaxDistanceThreshold << "\n"
			<< "--bfeature_validation " << bfeature_validation << "\n"
			<< "--bopticalfiltering " << bopticalfiltering << "\n"
			<< "--bdynamicdistance " << bdynamicdistance << "\n"
			<< "--bopticalmatching " << bopticalmatching << "\n"
			<< "--bSolveArticulationPoint " << bSolveArticulationPoint << "\n";
				
  EPairMode ePairmode = (iMatchingVideoMode == -1 ) ? PAIR_EXHAUSTIVE : PAIR_CONTIGUOUS;

  if (sPredefinedPairList.length()) {
    ePairmode = PAIR_FROM_FILE;
    if (iMatchingVideoMode>0) {
      std::cerr << "\nIncompatible options: --videoModeMatching and --pairList" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sMatchesDirectory.empty() || !stlplus::is_folder(sMatchesDirectory))  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }
  EGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
  std::string sGeometricMatchesFilename = "";
  switch (sGeometricModel[0])
  {
    case 'f': case 'F':
      eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
      sGeometricMatchesFilename = "matches.f.bin";
    break;
    case 'e': case 'E':
      eGeometricModelToCompute = ESSENTIAL_MATRIX;
      sGeometricMatchesFilename = "matches.e.bin";
    break;
    case 'h': case 'H':
      eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
      sGeometricMatchesFilename = "matches.h.bin";
    break;
    case 'a': case 'A':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ANGULAR;
      sGeometricMatchesFilename = "matches.f.bin";
    break;
    case 'o': case 'O':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ORTHO;
      sGeometricMatchesFilename = "matches.o.bin";
    break;
    default:
      std::cerr << "Unknown geometric model" << std::endl;
      return EXIT_FAILURE;
  }

  // -----------------------------
  // - Load SfM_Data Views & intrinsics data
  // a. Compute putative descriptor matches
  // b. Geometric filtering of putative matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // a. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider;
  if (ui_max_cache_size == 0)
  {
    // Default regions provider (load & store all regions in memory)
    regions_provider = std::make_shared<Regions_Provider>();
  }
  else
  {
    // Cached regions provider (load & store regions on demand)
    regions_provider = std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
  }

  // Show the progress on the command line:
  C_Progress_display progress;

  if (!regions_provider->load(sfm_data, sMatchesDirectory, regions_type, &progress)) {
    std::cerr << std::endl << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativesMatches;

  // Build some alias from SfM_Data Views data:
  // - List views as a vector of filenames & image sizes
  std::vector<std::string> vec_fileNames;
  std::vector<std::pair<size_t, size_t>> vec_imagesSize;
  {
    vec_fileNames.reserve(sfm_data.GetViews().size());
    vec_imagesSize.reserve(sfm_data.GetViews().size());
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
      iter != sfm_data.GetViews().end();
      ++iter)
    {
      const View * v = iter->second.get();
      vec_fileNames.push_back(stlplus::create_filespec(sfm_data.s_root_path,
          v->s_Img_path));
      vec_imagesSize.push_back( std::make_pair( v->ui_width, v->ui_height) );
    }
  }

  std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;
  // If the matches already exists, reload them
  if (!bForce
        && (stlplus::file_exists(sMatchesDirectory + "/matches.putative.txt")
        || stlplus::file_exists(sMatchesDirectory + "/matches.putative.bin"))
  )
  {
    if (!(Load(map_PutativesMatches, sMatchesDirectory + "/matches.putative.bin") ||
          Load(map_PutativesMatches, sMatchesDirectory + "/matches.putative.txt")) )
    {
      std::cerr << "Cannot load input matches file";
      return EXIT_FAILURE;
    }
    std::cout << "\t PREVIOUS RESULTS LOADED;"
      << " #pair: " << map_PutativesMatches.size() << std::endl;
  }
  else // Compute the putative matches
  {
    std::cout << "Use: ";
    switch (ePairmode)
    {
      case PAIR_EXHAUSTIVE: std::cout << "exhaustive pairwise matching" << std::endl; break;
      case PAIR_CONTIGUOUS: std::cout << "sequence pairwise matching" << std::endl; break;
      case PAIR_FROM_FILE:  std::cout << "user defined pairwise matching" << std::endl; break;
    }

    // Allocate the right Matcher according the Matching requested method
    std::unique_ptr<Matcher> collectionMatcher;
	std::unique_ptr<Hierarchical_Matcher_Regions> hierarchicalMatcher;
    if (sNearestMatchingMethod == "AUTO")
    {
      if (regions_type->IsScalar())
      {
        std::cout << "Using FAST_CASCADE_HASHING_L2 matcher" << std::endl;
        collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
      }
      else
      if (regions_type->IsBinary())
      {
        std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
        collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
      }
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEL2")
    {
      std::cout << "Using BRUTE_FORCE_L2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_L2));
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEHAMMING")
    {
      std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "HNSWL2")
    {
      std::cout << "Using HNSWL2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L2));
    }
    else
    if (sNearestMatchingMethod == "ANNL2")
    {
      std::cout << "Using ANN_L2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, ANN_L2));
    }
    else
    if (sNearestMatchingMethod == "CASCADEHASHINGL2")
    {
      std::cout << "Using CASCADE_HASHING_L2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, CASCADE_HASHING_L2));
    }
    else
    if (sNearestMatchingMethod == "FASTCASCADEHASHINGL2")
    {
      std::cout << "Using FAST_CASCADE_HASHING_L2 matcher" << std::endl;
      collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    }
	else if (sNearestMatchingMethod == "OPTICALFLOW" )
	{
		std::cout << "Using Optical_Flow_Matcher_Regions" << std::endl;
		collectionMatcher.reset(new Optical_Flow_Matcher_Regions(fDistRatio,bin_dir,MaxDistanceThreshold,sfm_data,sMatchesDirectory));
	}
	else if (sNearestMatchingMethod == "HIERARCHICAL")
	{
		std::cout << "Using Hierarchical_Matcher_Regions" << std::endl;
		hierarchicalMatcher.reset(
			new Hierarchical_Matcher_Regions(fDistRatio, bin_dir, MaxDistanceThreshold, sfm_data, sMatchesDirectory,
				bfeature_validation, bopticalfiltering, bdynamicdistance, bopticalmatching));
		//collectionMatcher.reset(new Hierarchical_Matcher_Regions(fDistRatio, bin_dir, MaxDistanceThreshold, sfm_data,sMatchesDirectory));
	}
    if (!collectionMatcher && !hierarchicalMatcher)
    {
      std::cerr << "Invalid Nearest Neighbor method: " << sNearestMatchingMethod << std::endl;
      return EXIT_FAILURE;
    }
    // Perform the matching
    system::Timer timer;
    {
      // From matching mode compute the pair list that have to be matched:
      Pair_Set pairs;
      switch (ePairmode)
      {
        case PAIR_EXHAUSTIVE: pairs = exhaustivePairs(sfm_data.GetViews().size()); break;
        case PAIR_CONTIGUOUS: pairs = contiguousWithOverlap(sfm_data.GetViews().size(), iMatchingVideoMode); break;
        case PAIR_FROM_FILE:
          if (!loadPairs(sfm_data.GetViews().size(), sPredefinedPairList, pairs))
          {
              return EXIT_FAILURE;
          }
          break;
      }
      // Photometric matching of putative pairs
	  if (sNearestMatchingMethod == "HIERARCHICAL")
	  {
		  hierarchicalMatcher->Hierarchical_Match(regions_provider, pairs, map_PutativesMatches, &progress);

		  //save the notable features
		  {
			  std::string output_dir(stlplus::create_filespec(sMatchesDirectory, "notables_features"));
			  stlplus::folder_create(output_dir);
			  std::ofstream notable_info(stlplus::create_filespec(output_dir, "notable_info.txt"));
			  for (Views::const_iterator iter = sfm_data.GetViews().begin();
				  iter != sfm_data.GetViews().end(); ++iter)
			  {
				  const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second->s_Img_path);
				  const std::string basename = stlplus::basename_part(sImageName);
				  const std::string featFile = stlplus::create_filespec(output_dir, basename, ".txt");
				  if (hierarchicalMatcher->notable_features.count(iter->first))
				  {
					  std::ofstream file(featFile);

					  std::shared_ptr<features::Regions> feats_1 = regions_provider->get(iter->first);
					  for (const auto& feature_id : hierarchicalMatcher->notable_features.at(iter->first))
					  {
						  Vec2 first_feat = feats_1->GetRegionPosition(feature_id);
						  file << feature_id << " " << first_feat(0) << " " << first_feat(1) << "\n";
					  }
					  file.close();
					  notable_info << basename << ":" << hierarchicalMatcher->notable_features.at(iter->first).size() << "\n";
				  }

			  }
			  notable_info.close();
		  }
	  }
	  else
	  {
		  collectionMatcher->Match(regions_provider, pairs, map_PutativesMatches, &progress);
	  }
      
      //---------------------------------------
      //-- Export putative matches
      //---------------------------------------
      if (!Save(map_PutativesMatches, std::string(sMatchesDirectory + "/matches.putative.bin")))
      {
        std::cerr
          << "Cannot save computed matches in: "
          << std::string(sMatchesDirectory + "/matches.putative.bin");
        return EXIT_FAILURE;
      }
    }
    std::cout << "Task (Regions Matching) done in (s): " << timer.elapsed() << std::endl;
  }
  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
    map_PutativesMatches,
    stlplus::create_filespec(sMatchesDirectory, "PutativeAdjacencyMatrix", "svg"));
  //-- export view pair graph once putative graph matches have been computed
  {
    std::set<IndexT> set_ViewIds;
    std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
      std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
    graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_PutativesMatches));
    graph::exportToGraphvizData(
      stlplus::create_filespec(sMatchesDirectory, "putative_matches"),
      putativeGraph);
  }

  //---------------------------------------
  // b. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------

  std::unique_ptr<ImageCollectionGeometricFilter> filter_ptr(
    new ImageCollectionGeometricFilter(&sfm_data, regions_provider));

  if (filter_ptr)
  {
    system::Timer timer;
    const double d_distance_ratio = 0.6;

    PairWiseMatches map_GeometricMatches;
    switch (eGeometricModelToCompute)
    {
      case HOMOGRAPHY_MATRIX:
      {
        const bool bGeometric_only_guided_matching = true;
        filter_ptr->Robust_model_estimation(
          GeometricFilter_HMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching,
          bGeometric_only_guided_matching ? -1.0 : d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case FUNDAMENTAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_FMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_EMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
        
        //-- Perform an additional check to remove pairs with poor overlap
        std::vector<PairWiseMatches::key_type> vec_toRemove;
        for (const auto & pairwisematches_it : map_GeometricMatches)
        {
          const size_t putativePhotometricCount = map_PutativesMatches.find(pairwisematches_it.first)->second.size();
          const size_t putativeGeometricCount = pairwisematches_it.second.size();
          const float ratio = putativeGeometricCount / static_cast<float>(putativePhotometricCount);
          if (putativeGeometricCount < 50 || ratio < .3f)  {
            // the pair will be removed
            vec_toRemove.push_back(pairwisematches_it.first);
          }
        }
        //-- remove discarded pairs
        for (const auto & pair_to_remove_it : vec_toRemove)
        {
          map_GeometricMatches.erase(pair_to_remove_it);
        }
      }
      break;
      case ESSENTIAL_MATRIX_ANGULAR:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_ESphericalMatrix_AC_Angular(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX_ORTHO:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_EOMatrix_RA(2.0, imax_iteration),
          map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
    }

	//////////Detect Articulation Point//////////////
	if (bSolveArticulationPoint)
	{

		//---------------------------------------
		// c. compute matches of pairs for eliminate the articulation point
		//---------------------------------------
		Pair_Set additional_pairs;
		PairWiseMatches art_PutativesMatches;
		SolveArticulationPoint(map_GeometricMatches, iMatchingVideoMode, additional_pairs);
		if (!additional_pairs.empty())
		{
			Cascade_Hashing_Matcher_Regions ch_matcher(fDistRatio);
			ch_matcher.Match(regions_provider, additional_pairs, art_PutativesMatches, &progress);
			/*Optical_Flow_Matcher_Regions of_matcher(fDistRatio, bin_dir, MaxDistanceThreshold, sfm_data, sMatchesDirectory);
			std::cout << "Solve Articulation Points\n";
			of_matcher.Match(regions_provider, additional_pairs, art_PutativesMatches, &progress);*/
			if (!art_PutativesMatches.empty())
			{
				for (const auto& pairwisematch : art_PutativesMatches)
				{
					if (!map_PutativesMatches.count(pairwisematch.first))
					{
						map_PutativesMatches.emplace(pairwisematch.first, pairwisematch.second);
					}
					else
					{
						std::vector<IndMatch>& indmatches = map_PutativesMatches.at(pairwisematch.first);
						std::set<IndexT> feat1_paired;
						std::set<IndexT> feat2_paired;
						for (const auto& indmatch_raw : indmatches)
						{
							feat1_paired.insert(indmatch_raw.i_);
							feat2_paired.insert(indmatch_raw.j_);
						}
						for (const auto& indmatch_extra : pairwisematch.second)
						{
							if (!feat1_paired.count(indmatch_extra.i_) && !feat2_paired.count(indmatch_extra.j_))
							{
								indmatches.emplace_back(indmatch_extra);
							}
						}
					}
				}

				//---------------------------------------
				// d. second Geometric filtering of  new matches
				//---------------------------------------

				std::unique_ptr<ImageCollectionGeometricFilter> filter_ptr(
					new ImageCollectionGeometricFilter(&sfm_data, regions_provider));

				if (filter_ptr)
				{


					map_GeometricMatches.clear();
					switch (eGeometricModelToCompute)
					{
					case HOMOGRAPHY_MATRIX:
					{
						const bool bGeometric_only_guided_matching = true;
						filter_ptr->Robust_model_estimation(
							GeometricFilter_HMatrix_AC(4.0, imax_iteration),
							map_PutativesMatches, bGuided_matching,
							bGeometric_only_guided_matching ? -1.0 : d_distance_ratio, &progress);
						map_GeometricMatches = filter_ptr->Get_geometric_matches();
					}
					break;
					case FUNDAMENTAL_MATRIX:
					{
						filter_ptr->Robust_model_estimation(
							GeometricFilter_FMatrix_AC(4.0, imax_iteration),
							map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
						map_GeometricMatches = filter_ptr->Get_geometric_matches();
					}
					break;
					case ESSENTIAL_MATRIX:
					{
						filter_ptr->Robust_model_estimation(
							GeometricFilter_EMatrix_AC(4.0, imax_iteration),
							map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
						map_GeometricMatches = filter_ptr->Get_geometric_matches();

						//-- Perform an additional check to remove pairs with poor overlap
						std::vector<PairWiseMatches::key_type> vec_toRemove;
						for (const auto & pairwisematches_it : map_GeometricMatches)
						{
							const size_t putativePhotometricCount = map_PutativesMatches.find(pairwisematches_it.first)->second.size();
							const size_t putativeGeometricCount = pairwisematches_it.second.size();
							const float ratio = putativeGeometricCount / static_cast<float>(putativePhotometricCount);
							if (putativeGeometricCount < 50 || ratio < .3f) {
								// the pair will be removed
								vec_toRemove.push_back(pairwisematches_it.first);
							}
						}
						//-- remove discarded pairs
						for (const auto & pair_to_remove_it : vec_toRemove)
						{
							map_GeometricMatches.erase(pair_to_remove_it);
						}
					}
					break;
					case ESSENTIAL_MATRIX_ANGULAR:
					{
						filter_ptr->Robust_model_estimation(
							GeometricFilter_ESphericalMatrix_AC_Angular(4.0, imax_iteration),
							map_PutativesMatches, bGuided_matching);
						map_GeometricMatches = filter_ptr->Get_geometric_matches();
					}
					break;
					case ESSENTIAL_MATRIX_ORTHO:
					{
						filter_ptr->Robust_model_estimation(
							GeometricFilter_EOMatrix_RA(2.0, imax_iteration),
							map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
						map_GeometricMatches = filter_ptr->Get_geometric_matches();
					}
					break;
					}
				}
			}
		}
	}
	

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    if (!Save(map_GeometricMatches,
      std::string(sMatchesDirectory + "/" + sGeometricMatchesFilename)))
    {
      std::cerr
          << "Cannot save computed matches in: "
          << std::string(sMatchesDirectory + "/" + sGeometricMatchesFilename);
      return EXIT_FAILURE;
    }

    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;

    // -- export Geometric View Graph statistics
    graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_GeometricMatches));

    //-- export Adjacency matrix
    std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
      << std::endl;
    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
      map_GeometricMatches,
      stlplus::create_filespec(sMatchesDirectory, "GeometricAdjacencyMatrix", "svg"));

    //-- export view pair graph once geometric filter have been done
    {
      std::set<IndexT> set_ViewIds;
      std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_GeometricMatches));
      graph::exportToGraphvizData(
        stlplus::create_filespec(sMatchesDirectory, "geometric_matches"),
        putativeGraph);
    }
  }
  return EXIT_SUCCESS;
}


