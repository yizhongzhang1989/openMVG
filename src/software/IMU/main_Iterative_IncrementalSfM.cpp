// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"
#include "openMVG/stl/stl.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

/// From 2 given image file-names, find the two corresponding index in the View list
bool computeIndexFromImageNames(
  const SfM_Data & sfm_data,
  const std::pair<std::string,std::string>& initialPairName,
  Pair& initialPairIndex)
{
  if (initialPairName.first == initialPairName.second)
  {
    std::cerr << "\nInvalid image names. You cannot use the same image to initialize a pair." << std::endl;
    return false;
  }

  initialPairIndex = {UndefinedIndexT, UndefinedIndexT};

  /// List views filenames and find the one that correspond to the user ones:
  for (Views::const_iterator it = sfm_data.GetViews().begin();
    it != sfm_data.GetViews().end(); ++it)
  {
    const View * v = it->second.get();
    const std::string filename = stlplus::filename_part(v->s_Img_path);
    if (filename == initialPairName.first)
    {
      initialPairIndex.first = v->id_view;
    }
    else{
      if (filename == initialPairName.second)
      {
        initialPairIndex.second = v->id_view;
      }
    }
  }
  return (initialPairIndex.first != UndefinedIndexT &&
      initialPairIndex.second != UndefinedIndexT);
}
////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
////////////////////
//find initial pairs associated with the first view in descending order of their correspondences count
//                          
//
// @param[in]  first_view                 the first view selected as initial pair
// @param[in]  map_Matches				  information of matches
// @param[in]  reindex_validview2index    the map from view id to discrete index
// @param[in]  tried_initialpairs         the initial pairs that have been tried once

// @param[out] initial_pairs              the initial pairs found

// (Owned by BC)
bool FindSecondImage(const IndexT first_view, const openMVG::matching::PairWiseMatches& map_Matches,
	const std::map<IndexT, IndexT>& reindex_validview2index,
	 const Pair_Set& tried_initialpairs,
	std::vector<Pair>& initial_pairs)
{
	//////////////////////////////////////////////
	///////////find second initial view///////////
	//////////////////////////////////////////////
	size_t second_image = 0;

	// list views matched with first_image
	std::vector<uint32_t> vec_NbMatchesPeView_J;
	std::vector<openMVG::matching::PairWiseMatches::const_iterator> vec_ViewIterator_J;
	for (openMVG::matching::PairWiseMatches::const_iterator
		iter = map_Matches.begin();
		iter != map_Matches.end(); ++iter)
	{
		const Pair current_pair = iter->first;
		if ((reindex_validview2index.count(current_pair.first) &&
			reindex_validview2index.count(current_pair.second)) &&
			!tried_initialpairs.count(current_pair) &&
			(current_pair.first == first_view || current_pair.second == first_view))
		{
			vec_NbMatchesPeView_J.push_back(iter->second.size());
			vec_ViewIterator_J.push_back(iter);
		}
	}

	// sort the Views in descending order of correspondences
	using namespace stl::indexed_sort;
	std::vector<sort_index_packet_descend<uint32_t, uint32_t>> packet_vec_J(vec_NbMatchesPeView_J.size());
	sort_index_helper(packet_vec_J, &vec_NbMatchesPeView_J[0]);


	//return the corresponds with first_view in descending order of their correspondences count
	for (size_t i = 0; i < vec_NbMatchesPeView_J.size(); ++i)
	{
		const uint32_t index = packet_vec_J[i].index;
		openMVG::matching::PairWiseMatches::const_iterator iter = vec_ViewIterator_J[index];
		initial_pairs.push_back(iter->first);
	}

	return initial_pairs.size() > 0;
}
////////////////////
//find the initial pair according the corresponds
//    * the first view is selected as the largest connected view among unregistered views
//    * the second views are the views connected with first view ordered in descending of correspondences.
//                          
//
// @param[in]  sfm_data                   input sfm data to provide views' information
// @param[in]  registered_views           the view pairs to matching
// @param[in]  matches_provider           matching information
// @param[in]  tried_initialpairs         the initial pairs that have been tried once
// @param[in]  first_initial_tried_images the views that have been selected as first view of initial pair
// @param[out] initial_pairs              the intial pairs found

// (Owned by BC)
bool FindInitialImagePairs(const SfM_Data& sfm_data,const std::set<IndexT>& registered_views,
                           const Matches_Provider* matches_provider,const Pair_Set& tried_initialpairs,
						  std::set<IndexT>& first_initial_tried_images,
						   std::vector<Pair>& initial_pairs)
{

    
    // List Views that supports valid intrinsic and are not yet registered
    std::vector<uint32_t> vec_NbMatchesPerView_I;
    std::vector<IndexT> vec_ViewsIterator_I;
    vec_NbMatchesPerView_I.reserve(sfm_data.GetViews().size());
    vec_ViewsIterator_I.reserve(sfm_data.GetViews().size());

	std::map<IndexT, IndexT> reindex_validview2index;  	//used for discretization
	IndexT reindex_cpt = 0;
    for (const auto & view : sfm_data.GetViews())
    {
      const View * v = view.second.get();
      if (sfm_data.GetIntrinsics().find(v->id_intrinsic) != sfm_data.GetIntrinsics().end()
          && !registered_views.count(view.first))
      {
        reindex_validview2index.emplace(v->id_view,reindex_cpt);
        vec_NbMatchesPerView_I.push_back(0);
        vec_ViewsIterator_I.push_back(v->id_view);
        reindex_cpt++;
        
      }
    }
    ////////////////////////////////////////////////////
    ///////////select first initial view////////////////
    ////////////////////////////////////////////////////
    IndexT first_view = 0;
    
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider->pairWise_matches_;

    // count the number of correspondings per views
    for (openMVG::matching::PairWiseMatches::const_iterator
      iter = map_Matches.begin(); iter != map_Matches.end(); ++iter)
    {
      const Pair current_pair = iter->first;
      if (reindex_validview2index.count(current_pair.first) &&
        reindex_validview2index.count(current_pair.second) )
      {
        const IndexT& index_1 = reindex_validview2index.at(current_pair.first);
        const IndexT& index_2 = reindex_validview2index.at(current_pair.second);
        vec_NbMatchesPerView_I[index_1] += iter->second.size();
        vec_NbMatchesPerView_I[index_2] += iter->second.size();
      }
    }

    // sort the Views in descending order of their correspondences count
    using namespace stl::indexed_sort;
    std::vector<sort_index_packet_descend<uint32_t, uint32_t>> packet_vec_I(vec_NbMatchesPerView_I.size());
    sort_index_helper(packet_vec_I, &vec_NbMatchesPerView_I[0]);


    //select the largest connected view as first view
    bool find_firstview = false;
    for (size_t i = 0; i < vec_NbMatchesPerView_I.size(); ++i) 
    {
        const uint32_t index = packet_vec_I[i].index;
        IndexT view_i = vec_ViewsIterator_I[index];
		initial_pairs.clear();
		assert(initial_pairs.size() == 0);
        if(!first_initial_tried_images.count(view_i)  // the view i has not been selected as first initial image.
			&& FindSecondImage(view_i,map_Matches,
				reindex_validview2index,
				tried_initialpairs,
				initial_pairs)) {
          first_view = view_i;
          find_firstview = true;
          break;
        }
    }
    if(!find_firstview)
    {
        return false;
    }

    first_initial_tried_images.insert(first_view);

	return true;
    


}
//END(Author: BC)===================================================


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Iterative Sequential/Incremental reconstruction" << std::endl
            << " Perform incremental SfM (Initial Pair Essential + Resection)." << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir, sMatchFilename;
  std::string sOutDir = "";
  std::pair<std::string,std::string> initialPairString("","");
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
  bool b_use_motion_priors = false;
  int triangulation_method = static_cast<int>(ETriangulationMethod::DEFAULT);

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('M', sMatchFilename, "match_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('a', initialPairString.first, "initialPairA") );
  cmd.add( make_option('b', initialPairString.second, "initialPairB") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_option('t', triangulation_method, "triangulation_method"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n[Optional]\n"
    << "[-a|--initialPairA] filename of the first image (without path)\n"
    << "[-b|--initialPairB] filename of the second image (without path)\n"
    << "[-c|--camera_model] Camera model type for view with unknown intrinsic:\n"
      << "\t 1: Pinhole \n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole radial 3 + tangential 2\n"
      << "\t 5: Pinhole fisheye\n"
    << "[-f|--refineIntrinsics] Intrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> intrinsic parameters are held as constant\n"
      << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
      << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
      << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
      << "\t -> NOTE: options can be combined thanks to '|'\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
      <<      "\t\t-> refine the focal length & the principal point position\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
      << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
    << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions) (default: false)\n"
    << "[-M|--match_file] path to the match file to use (default=matches.f.txt then matches.f.bin).\n"
    << "[-t|--triangulation_method] triangulation method (default=" << triangulation_method << "):\n"
    << "\t\t" << static_cast<int>(ETriangulationMethod::DIRECT_LINEAR_TRANSFORM) << ": DIRECT_LINEAR_TRANSFORM\n"
    << "\t\t" << static_cast<int>(ETriangulationMethod::L1_ANGULAR) << ": L1_ANGULAR\n"
    << "\t\t" << static_cast<int>(ETriangulationMethod::LINFINITY_ANGULAR) << ": LINFINITY_ANGULAR\n"
    << "\t\t" << static_cast<int>(ETriangulationMethod::INVERSE_DEPTH_WEIGHTED_MIDPOINT) << ": INVERSE_DEPTH_WEIGHTED_MIDPOINT\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if ( !isValid(static_cast<ETriangulationMethod>(triangulation_method))) {
    std::cerr << "\n Invalid triangulation method" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !isValid(openMVG::cameras::EINTRINSIC(i_User_camera_model)) )  {
    std::cerr << "\n Invalid camera type" << std::endl;
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
    cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
  {
    std::cerr << "Invalid input for Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data raw_sfm_data;
  if (!Load(raw_sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(raw_sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  // Matches reading
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if // Try to read the provided match filename or the default one (matches.f.txt/bin)
  (
    !(matches_provider->load(raw_sfm_data, sMatchFilename) ||
      matches_provider->load(raw_sfm_data, stlplus::create_filespec(sMatchesDir, "matches.f.txt")) ||
      matches_provider->load(raw_sfm_data, stlplus::create_filespec(sMatchesDir, "matches.f.bin")))
  )
  {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "\nCannot create the output directory" << std::endl;
    }
  }

  //---------------------------------------
  // Iterative Sequential reconstruction process
  //---------------------------------------
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  openMVG::system::Timer timer;
  std::set<IndexT> registered_views;            //record the views have been registered

  std::set<IndexT> first_initial_tried_images;  //record the views have been tried as first view
  Pair_Set tried_initialpairs;                  //record the initial pairs have been tried once
  int iteration_i = 1;                          
                   
  const size_t num_views = raw_sfm_data.GetViews().size();
  size_t num_reconstructions = 0;
  while(num_views - registered_views.size() > 3)
  {
    std::cout<<"/////////////////////////\n"
             <<"///////iteration "<<iteration_i<<"///////\n"
             <<"/////////////////////////\n";
    
    // find the initial pairs among  unregistered views
    std::vector<Pair> initial_pairs;
    if(iteration_i==1)  //In the first iteration ,we let OpenMVG automatically find initial pairs.
    {
        initial_pairs.push_back(Pair(0,0)); 
    }
    else
    {
      if(!FindInitialImagePairs(raw_sfm_data,registered_views,
                                matches_provider.get(),tried_initialpairs, 
                                first_initial_tried_images,
                                initial_pairs))
      {
        std::cout<<"There are not good initial pairs for next reconstruction\n";
        break;
      }
    }

	//create output directory
	std::string iteration_sOutDir = stlplus::create_filespec(sOutDir, std::to_string(iteration_i));
	if (!stlplus::folder_exists(iteration_sOutDir))
	{
		if (!stlplus::folder_create(iteration_sOutDir))
		{
			std::cerr << "\nCannot create the output directory" << std::endl;
		}
	}

    //try to reconstruction with initial pair
    //  * If reconstruct successfully,exit loop and find next initial pairs for unregistered views
	bool brecons = false;
    for(const Pair& initial_pair:initial_pairs)
    {
      openMVG::system::Timer re_timer;
      SequentialSfMReconstructionEngine sfmEngine(
        raw_sfm_data,                             //pass by value
        iteration_sOutDir,
        stlplus::create_filespec(iteration_sOutDir, "Reconstruction_Report.html"));

      // Configure the features_provider & the matches_provider
      sfmEngine.SetFeaturesProvider(feats_provider.get());
      sfmEngine.SetMatchesProvider(matches_provider.get());

      // Configure reconstruction parameters
      sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
      sfmEngine.SetUnknownCameraType(EINTRINSIC(i_User_camera_model));
      b_use_motion_priors = cmd.used('P');
      sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);
      sfmEngine.SetTriangulationMethod(static_cast<ETriangulationMethod>(triangulation_method));

      if(iteration_i == 1 && initial_pair.first==0&& initial_pair.second == 0)
      {
        //User can specify the initial pair in the first reconstruction
        // Handle Initial pair parameter
        if (!initialPairString.first.empty() && !initialPairString.second.empty())
        {
          Pair initialPairIndex;
          if (!computeIndexFromImageNames(raw_sfm_data, initialPairString, initialPairIndex))
          {
              std::cerr << "Could not find the initial pairs <" << initialPairString.first
                <<  ", " << initialPairString.second << ">!\n";
            return EXIT_FAILURE;
          }
          sfmEngine.setInitialPair(initialPairIndex);
        }
      }
      else
      {
		 std::cout << "reconstruction begins with view pair (" << initial_pair.first << "," << initial_pair.second << ")\n";
         sfmEngine.setInitialPair(initial_pair);
         tried_initialpairs.insert(initial_pair);   //record the initial pairs used.
      }
	  //start reconstruction
	  bool bProcess = sfmEngine.Process();
	   
      if (bProcess && sfmEngine.Get_SfM_Data().GetPoses().size() > 3)
      {
		const SfM_Data& iteration_sfm_data = sfmEngine.Get_SfM_Data();
		brecons = true;
        num_reconstructions ++;
        std::cout << "...Generating SfM_Report.html" << std::endl;
        Generate_SfM_Report(iteration_sfm_data,
          stlplus::create_filespec(iteration_sOutDir, "SfMReconstruction_Report.html"));

        //-- Export to disk computed scene (data & visualizable results)
        std::cout << "...Export SfM_Data to disk." << std::endl;
        Save(iteration_sfm_data,
          stlplus::create_filespec(iteration_sOutDir, "sfm_data", ".bin"),
          ESfM_Data(ALL));

        Save(iteration_sfm_data,
          stlplus::create_filespec(iteration_sOutDir, "cloud_and_poses", ".ply"),
          ESfM_Data(ALL));

		//assure the value of raw_sfm_data is not modified.
        if(raw_sfm_data.GetPoses().size() > 0)
        {
          std::cout<<"Error:raw_sfm_data has been modified\n";
          return EXIT_SUCCESS;
        }

        //insert the views registered into `registered_views`
        
        size_t newviews = 0;
        for(const auto& view_item : iteration_sfm_data.GetViews())
        {
            if(iteration_sfm_data.IsPoseAndIntrinsicDefined(view_item.second.get()))
            {
               if (!registered_views.count(view_item.first)) newviews ++;
               registered_views.insert(view_item.first);
            }
        }
        
        std::cout<<"Static/////////////////////\n";
        std::cout<<"#Images in this reconstruction:"<<iteration_sfm_data.GetPoses().size()<<"\n";
        std::cout<<"#Images newly registered:"<<newviews<<"\n";
        std::cout<<"#Images registered:"<<registered_views.size()<<"\n";
        std::cout<<"#Images not yet registered:"<<num_views - registered_views.size()<<"\n";
        std::cout << std::endl << "Ac-Sfm took (s): " << re_timer.elapsed() << std::endl;
        std::cout<<"/////////////////////\n";
        break;
        
      }

    }
	//delete the directory already created.
	if (!brecons)
	{
		if (!stlplus::folder_delete(iteration_sOutDir, true))
		{
			std::cout << "Error:Can not delete directory:" << iteration_sOutDir << "\n";
			return EXIT_FAILURE;
		}
	}
	else
	{
		iteration_i++;
	}
  }
  std::cout<<"\n////////////Holistic Static//////////////\n";
  std::cout<<"#Images registered:"<<registered_views.size()<<"\n";
  std::cout<<"#Images not yet registered:"<<num_views - registered_views.size()<<"\n";
  std::cout << std::endl << " Total Ac-Sfm took (s): " <<timer.elapsed() << std::endl;
  std::cout<<"/////////////////////\n";
  //END(Author: BC)===================================================
  return EXIT_SUCCESS;
}

