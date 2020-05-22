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

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

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

inline int double_comp(double x, double y)
{
	double eps = 1e-5;
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

inline double cal_dis(double x1, double y1, double x2, double y2)
{
	return std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}


bool OpticalFiltering(PairWiseMatches& map_PutativesMatches,const Regions_Provider* regions_provider,
					  const std::string bin_dir,const SfM_Data& sfm_data,const std::string matches_dir
					  ,const double MaxDistanceThreshold)
{
  // std::shared_ptr<features::Regions> feat_0= regions_provider->get(0);
  // const features::PointFeatures feats = feat_0->GetRegionsPositions();
  // for(size_t i = 0; i<feats.size(); i++)
  // {
  //     std::cout<<"feat "<<i<<" in view 0:("<<feats[i].x()<<","<<feats[i].y()<<") \n";
  // }
  std::cout << "*****Optical Filtering****\n";
  struct Optical_track
  {
    unsigned int view_id;
		unsigned int feat_id;
		float x;
		float y;
		float rx;
		float ry;
  };
  std::map<IndexT,std::map<Pair,Optical_track>> opticaltrack_table;
  //read binary data of optical flow for every features
  std::cout << "**read binary features**\n";
  for (const auto& view_item:sfm_data.GetViews())
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
		  return false;
	  }
      if(opticaltrack_table.count(view_id))
      {
        std::cerr<<"Error:duplicate view id\n";
        return false;
      }
      opticaltrack_table.emplace(view_id,std::map<Pair,Optical_track>());
      std::map<Pair,Optical_track>& map_optrack = opticaltrack_table.at(view_id);

      while (file.read(reinterpret_cast<char*>(&ot), sizeof(Optical_track)))
      {
        map_optrack[Pair(ot.view_id,ot.feat_id)] = ot;
        if (file.peek() == '\n')
        {
          file.ignore();
        }
      }
	  std::cout << "length of view " << view_id << ": " << opticaltrack_table.at(view_id).size() << "\n";
      
  }
  //debugs//
  bool debug = false;
  IndexT img1 = 31;
  IndexT img2 = 36;
  std::ofstream filelog(stlplus::create_filespec(matches_dir, "debug_"+std::to_string(img1)+"_"+ std::to_string(img2) +".txt"));
  std::vector<double > dissss_;
  //debuge//
  //filter putative matches
  std::cout << "**filter putative matches**\n";
  PairWiseMatches::iterator pwm_iter;
  std::vector<double> dis_metric;
  size_t count_featpairs_untracked = 0;
  size_t count_featpairs_accepted = 0;
  size_t count_featpairs_all = 0;
  size_t count_featpairs_filtered = 0;
  for(pwm_iter=map_PutativesMatches.begin(); pwm_iter!=map_PutativesMatches.end(); pwm_iter++)
  {
      const IndexT first_view_id = pwm_iter->first.first;
      const IndexT second_view_id = pwm_iter->first.second;
	  const IndexT latter_view_id = std::max(first_view_id, second_view_id);
      std::shared_ptr<features::Regions> feats_1= regions_provider->get(first_view_id);
      std::shared_ptr<features::Regions> feats_2= regions_provider->get(second_view_id);
	  
	  std::cout << "Pair " << first_view_id << "-" << second_view_id << "\n";
	  if ((first_view_id == img1 && second_view_id == img2))
	  {
		  debug = true;
	  }
	  else
	  {
		  debug = false;
	  }
	  IndMatches::iterator im_iter;
	  const std::map<Pair, Optical_track>& latter_opticaltable = opticaltrack_table.at(latter_view_id);
	  std::cout << "	#raw feature pairs" << pwm_iter->second.size() << "\n";
	  for (im_iter = pwm_iter->second.begin(); im_iter!=pwm_iter->second.end();)
	  {
		  count_featpairs_all++;
		  const IndexT first_feat_id = im_iter->i_;
		  const IndexT second_feat_id = im_iter->j_;
		  Vec2 first_feat = feats_1->GetRegionPosition(first_feat_id);
		  Vec2 second_feat = feats_2->GetRegionPosition(second_feat_id);
		  if (latter_opticaltable.count(Pair(first_view_id, first_feat_id)) == 0 )
		  {
			  if (latter_view_id == first_view_id)
			  {
				  //the first feature  is not tracked in second view
				  std::cerr << "Error: feature" << first_view_id << " is losed in view " << first_view_id << "\n";
				  return false;
			  }
			  else
			  {
				 /* if (debug)
				  {
					  filelog << "feature " << first_feat_id << "-" << second_feat_id << " is not be tracked\n";
				  }*/
				  //the first feature  is not tracked in second view
				  im_iter = pwm_iter->second.erase(im_iter);
				  count_featpairs_untracked++;
				  continue;
			  }
			  
		  }
		  if (latter_opticaltable.count(Pair(second_view_id, second_feat_id)) == 0 )
		  {
			  if (latter_view_id == second_view_id)
			  {
				  //the first feature  is not tracked in second view
				  std::cerr << "Error: feature" << second_feat_id << " is losed in view " << second_view_id << "\n";
				  return false;
			  }
			  else
			  {
				  /*if (debug)
				  {
					  filelog << "feature " << first_feat_id << "-" << second_feat_id << " is not be tracked\n";
				  }*/
				  //the first feature  is not tracked in second view
				  im_iter = pwm_iter->second.erase(im_iter);
				  count_featpairs_untracked++;
				  continue;
			  }
		  }
			//check the distance of features in second view
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
			
			double dis_featpair = cal_dis(ot_1.x, ot_1.y, ot_2.x, ot_2.y);
			dis_metric.push_back(dis_featpair);
			if (dis_featpair > MaxDistanceThreshold)   //delete the feature pairs whose distance is too large.
			{
				if (debug)
				{
					filelog << ot_1.rx<<" "<<ot_1.ry<<" "<<ot_1.x << " " << ot_1.y << " " 
						<< ot_2.rx << " " << ot_2.ry << " " << ot_2.x << " " << ot_2.y << " "<< dis_featpair <<" 0\n";
					dissss_.push_back(dis_featpair);
				}
				count_featpairs_filtered++;
				im_iter = pwm_iter->second.erase(im_iter);
			}
			else
			{
				if (debug)
				{
					filelog << ot_1.rx << " " << ot_1.ry << " " << ot_1.x << " " << ot_1.y << " "
						<< ot_2.rx << " " << ot_2.ry << " " << ot_2.x << " " << ot_2.y << " " << dis_featpair << " 1\n";
					dissss_.push_back(dis_featpair);
				}
				count_featpairs_accepted++;
				im_iter++;
			}
			
			
      }
	  std::cout << "	#feature pairs accepted:" << pwm_iter->second.size() << "\n";
      
  }
  std::cout << "# feature pairs untracked:" << count_featpairs_untracked << "\n";
  std::cout << "# feature pairs accepted:" << count_featpairs_accepted << "\n";
  std::cout << "# feature pairs filtered:" << count_featpairs_filtered << "\n";
  std::cout << "# feature pairs all:" << count_featpairs_all << "\n";

  //debugs//
  if (dissss_.size()!=0)
  {
	  std::sort(dissss_.begin(), dissss_.end());
	  std::ofstream dissss_stream(stlplus::create_filespec(matches_dir, "dissss_"+std::to_string(img1)+"_"+ std::to_string(img2) +".txt"));
	  double sum = 0.0;
	  for (size_t i = 0; i < dissss_.size(); i++)
	  {
		  dissss_stream << dissss_[i] << "\n";
		  sum += dissss_[i];
	  }

	  dissss_stream << "Min:" << dissss_[0] << "\n";
	  dissss_stream << "Max:" << *(dissss_.rbegin()) << "\n";
	  dissss_stream << "Median:" << dissss_[dissss_.size() / 2] << "\n";
	  dissss_stream << "Mean:" << sum / ((double)dissss_.size()) << "\n";

	  dissss_stream << "# feature pairs untracked:" << count_featpairs_untracked << "\n";
	  dissss_stream << "# feature pairs accepted:" << count_featpairs_accepted << "\n";
	  dissss_stream << "# feature pairs filtered:" << count_featpairs_filtered << "\n";
	  dissss_stream << "# feature pairs all:" << count_featpairs_all << "\n";
	  dissss_stream.close();
  }
  filelog.close();
  //debuge//
  //output distances
  std::cout << "**output distances**\n";
  std::sort(dis_metric.begin(), dis_metric.end());
  std::ofstream dis_metric_stream(stlplus::create_filespec(matches_dir, "dis_metric.txt"));
  double sum = 0.0;
  for (size_t i = 0; i < dis_metric.size(); i++)
  {
	  dis_metric_stream << dis_metric[i] << "\n";
	  sum += dis_metric[i];
  }
  
  dis_metric_stream << "Min:" << dis_metric[0] << "\n";
  dis_metric_stream << "Max:" << *(dis_metric.rbegin()) << "\n";
  dis_metric_stream << "Median:" << dis_metric[dis_metric.size() / 2] << "\n";
  dis_metric_stream << "Mean:" << sum / ((double)dis_metric.size()) << "\n";

  dis_metric_stream << "# feature pairs untracked:" << count_featpairs_untracked << "\n";
  dis_metric_stream << "# feature pairs accepted:" << count_featpairs_accepted << "\n";
  dis_metric_stream << "# feature pairs filtered:" << count_featpairs_filtered << "\n";
  dis_metric_stream << "# feature pairs all:" << count_featpairs_all << "\n";
  dis_metric_stream.close();
}

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
/// - Export computed data
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
  double MaxDistanceThreshold = 40.0;

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
      << "  For Binary based descriptor:\n"
      << "    BRUTEFORCEHAMMING: BruteForce Hamming matching.\n"
      << "[-m|--guided_matching]\n"
      << "  use the found model to improve the pairwise correspondences.\n"
      << "[-c|--cache_size]\n"
      << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
      << "  If not used, all regions will be load in memory."
      << "[-b|--bin_dir]\n"
	    << "  the directory stores the binary file of position of features tracked in optical flow.\n"
	    << "  It is required for the optical filtering.\n"
      << "[-k|--maxdistancethreshold]\n"
	    << "  the distance threshold for the optical filtering.\n"
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
            << "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << std::endl;

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
    if (!collectionMatcher)
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
      collectionMatcher->Match(regions_provider, pairs, map_PutativesMatches, &progress);
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
  // b. Optical filtering of putative matches
  //    
  //    
  //---------------------------------------
  OpticalFiltering(map_PutativesMatches,regions_provider.get(), 
				   bin_dir, sfm_data,sMatchesDirectory, MaxDistanceThreshold);
  //---------------------------------------
  //-- Export putative matches
  //---------------------------------------
  if (!Save(map_PutativesMatches, std::string(sMatchesDirectory + "/matches.opticalfiltering.bin")))
  {
    std::cerr
      << "Cannot save computed matches in: "
      << std::string(sMatchesDirectory + "/matches.opticalfiltering.bin");
    return EXIT_FAILURE;
  }
  {
    std::set<IndexT> set_ViewIds;
    std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
      std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
    graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_PutativesMatches));
    graph::exportToGraphvizData(
      stlplus::create_filespec(sMatchesDirectory, "opticalfiltering"),
      putativeGraph);
  }
  //return EXIT_SUCCESS;
  //---------------------------------------
  // c. Geometric filtering of putative matches
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


