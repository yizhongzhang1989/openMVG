// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"

#include "openMVG_IMU/utils/myoutput.hpp"
#include "openMVG_IMU/sfm/pipelines/sfm_robust_multimodel_estimation.hpp"

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

std::shared_ptr<openMVG::tracks::SharedTrackVisibilityHelper> InitLandmarkTracks(const Matches_Provider* matches_provider, openMVG::tracks::STLMAPTracks& map_tracks_)
{
	// Compute tracks from matches
	tracks::TracksBuilder tracksBuilder;
	
	{
		// List of features matches for each couple of images
		const openMVG::matching::PairWiseMatches & map_Matches = matches_provider->pairWise_matches_;
		std::cout << "\n" << "Track building" << std::endl;

		tracksBuilder.Build(map_Matches);
		std::cout << "\n" << "Track filtering" << std::endl;
		tracksBuilder.Filter();
		std::cout << "\n" << "Track export to internal struct" << std::endl;
		//-- Build tracks with STL compliant type :
		
		tracksBuilder.ExportToSTL(map_tracks_);

		std::cout << "\n" << "Track stats" << std::endl;
		{
			std::ostringstream osTrack;
			//-- Display stats :
			//    - number of images
			//    - number of tracks
			std::set<uint32_t> set_imagesId;
			tracks::TracksUtilsMap::ImageIdInTracks(map_tracks_, set_imagesId);
			osTrack << "------------------" << "\n"
				<< "-- Tracks Stats --" << "\n"
				<< " Tracks number: " << tracksBuilder.NbTracks() << "\n"
				<< " Images Id: " << "\n";
			std::copy(set_imagesId.begin(),
				set_imagesId.end(),
				std::ostream_iterator<uint32_t>(osTrack, ", "));
			osTrack << "\n------------------" << "\n";

			std::map<uint32_t, uint32_t> map_Occurence_TrackLength;
			tracks::TracksUtilsMap::TracksLength(map_tracks_, map_Occurence_TrackLength);
			osTrack << "TrackLength, Occurrence" << "\n";
			for (const auto & it : map_Occurence_TrackLength) {
				osTrack << "\t" << it.first << "\t" << it.second << "\n";
			}
			osTrack << "\n";
			std::cout << osTrack.str();
		}
	}
	// Initialize the shared track visibility helper
	//shared_track_visibility_helper.reset(new openMVG::tracks::SharedTrackVisibilityHelper(map_tracks_));
	std::shared_ptr<openMVG::tracks::SharedTrackVisibilityHelper> shared_track_visibility_helper = std::make_shared<openMVG::tracks::SharedTrackVisibilityHelper>(map_tracks_);
	return shared_track_visibility_helper;
}


bool FindInitialImagePairs(const SfM_Data& sfm_data,const std::set<IndexT>& registered_views,
                           const Matches_Provider* matches_provider,Features_Provider* features_provider,const Pair_Set& tried_initialpairs,
						   Hash_Map<Pair,Pose3>& relative_poses, std::vector<std::pair<double, Pair>>& scoring_per_pair,
                           const std::string sOutDir, const int& iteration_i,const std::string file_suffix,const size_t max_iteration_count,
						   openMVG::tracks::SharedTrackVisibilityHelper* shared_track_visibility_helper)
{

	// select a pair that have the largest baseline (mean angle between its bearing vectors).

	const unsigned iMin_inliers_count = 100;
	const float fRequired_min_angle = 3.0f;
	const float fLimit_max_angle = 60.0f; // More than 60 degree, we cannot rely on matches for initial pair seeding

	
	
	// List Views that support valid intrinsic (view that could be used for Essential matrix computation)
	std::set<IndexT> valid_views;
	for (Views::const_iterator it = sfm_data.GetViews().begin();
		it != sfm_data.GetViews().end(); ++it)
	{
		const View * v = it->second.get();
		if (sfm_data.GetIntrinsics().count(v->id_intrinsic))
			valid_views.insert(v->id_view);
	}

	if (valid_views.size() < 2)
	{
		return false; // There is not view that support valid intrinsic data
	}
	//////////BC DEBUG START///////
	bool bdebug = true;
	std::ofstream essentialmatrix_log;
	std::ofstream homographymatrix_log;
	if(bdebug)
	{
		essentialmatrix_log.open(stlplus::create_filespec(sOutDir,"total_initial_essential_matrix.txt"));
		homographymatrix_log.open(stlplus::create_filespec(sOutDir, "total_initial_homography_matrix.txt"));
	}
	//////////BC DEBUG END///////
	// Compute the relative pose & the 'baseline score'
	C_Progress_display my_progress_bar(matches_provider->pairWise_matches_.size(),
		std::cout,
		"MY automatic selection of an initial pair:\n");
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
	for (const std::pair<Pair, matching::IndMatches> & match_pair : matches_provider->pairWise_matches_)
	{
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
		{
			++my_progress_bar;
			
			const Pair current_pair = match_pair.first;

			const uint32_t I = std::min(current_pair.first, current_pair.second);
			const uint32_t J = std::max(current_pair.first, current_pair.second);
			if (valid_views.count(I) && valid_views.count(J))
			{
				const View
					* view_I = sfm_data.GetViews().at(I).get(),
					*view_J = sfm_data.GetViews().at(J).get();
				const Intrinsics::const_iterator
					iterIntrinsic_I = sfm_data.GetIntrinsics().find(view_I->id_intrinsic),
					iterIntrinsic_J = sfm_data.GetIntrinsics().find(view_J->id_intrinsic);

				const auto
					cam_I = iterIntrinsic_I->second.get(),
					cam_J = iterIntrinsic_J->second.get();
				if (cam_I && cam_J)
				{
					openMVG::tracks::STLMAPTracks map_tracksCommon;
					shared_track_visibility_helper->GetTracksInImages({ I, J }, map_tracksCommon);
					// Copy points correspondences to arrays for relative pose estimation
					const size_t n = map_tracksCommon.size();
					Mat xI(2, n), xJ(2, n);
					size_t cptIndex = 0;
					for (const auto & track_iter : map_tracksCommon)
					{
						auto iter = track_iter.second.cbegin();
						const uint32_t i = iter->second;
						const uint32_t j = (++iter)->second;

						Vec2 feat = features_provider->feats_per_view[I][i].coords().cast<double>();
						xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
						feat = features_provider->feats_per_view[J][j].coords().cast<double>();
						xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
						++cptIndex;
					}

					// Robust estimation of the relative pose
					RelativePose_MultiInfo relativePose_info;
					RelativePose_MultiInfo relativePose_essentialinfo;
					RelativePose_MultiInfo relativePose_homographyinfo;
					//relativePose_info.initial_residual_tolerance = Square(4.0);   //(useless)can not send into function
					//relativePose_essentialinfo.initial_residual_tolerance = Square(4.0);  
					//relativePose_homographyinfo.initial_residual_tolerance = Square(4.0);

					if (robustRelativePose_MultiModel(
						cam_I, cam_J,
						xI, xJ, relativePose_info,
						{ cam_I->w(), cam_I->h() }, { cam_J->w(), cam_J->h() },
						max_iteration_count)
						&& relativePose_info.vec_inliers.size() > iMin_inliers_count)
					{
						// Triangulate inliers & compute angle between bearing vectors
						std::vector<float> vec_angles;
						vec_angles.reserve(relativePose_info.vec_inliers.size());
						const Pose3 pose_I = Pose3(Mat3::Identity(), Vec3::Zero());
						const Pose3 pose_J = relativePose_info.relativePose;
						
						for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
						{
							openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
							std::advance(iterT, inlier_idx);
							tracks::submapTrack::const_iterator iter = iterT->second.begin();
							const Vec2 featI = features_provider->feats_per_view[I][iter->second].coords().cast<double>();
							const Vec2 featJ = features_provider->feats_per_view[J][(++iter)->second].coords().cast<double>();
							vec_angles.push_back(AngleBetweenRay(pose_I, cam_I, pose_J, cam_J,
								cam_I->get_ud_pixel(featI), cam_J->get_ud_pixel(featJ)));
						}
						// Compute the median triangulation angle
						const unsigned median_index = vec_angles.size() / 2;
						std::nth_element(
							vec_angles.begin(),
							vec_angles.begin() + median_index,
							vec_angles.end());
						const float scoring_angle = vec_angles[median_index];
						// Store the pair iff the pair is in the asked angle range [fRequired_min_angle;fLimit_max_angle]
						if (scoring_angle > fRequired_min_angle &&
							scoring_angle < fLimit_max_angle)
						{
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
							{
								if(bdebug)  //bc
								{
									essentialmatrix_log << "Pair "<<current_pair.first<<" - "<< current_pair.second<<"\n";
									essentialmatrix_log <<"relation: "<<(relativePose_info.b_coplanar?"homography":"essential")<<"\n";
									
									
									essentialmatrix_log <<"essential:"<< relativePose_essentialinfo.vec_inliers.size()<<" in "<< xI .cols()<<"\n";
									
									for(int i = 0 ;i < 3; i++)
									{
										for(int j = 0 ;j < 3; j++)
										{
											essentialmatrix_log << relativePose_essentialinfo.model_matrix(i,j)<<" ";
										}
										essentialmatrix_log <<"\n";
									}
									essentialmatrix_log << "relative pose:\n";
									Mat3 R_E = relativePose_essentialinfo.relativePose.rotation();
									Vec3 t_E = relativePose_essentialinfo.relativePose.translation();
									for (int i = 0; i < 3; i++)
									{
										for (int j = 0; j < 3; j++)
										{
											essentialmatrix_log << R_E(i, j) << " ";
										}
										essentialmatrix_log << t_E(i) << "\n";
									}

									homographymatrix_log << "Pair " << current_pair.first << " - " << current_pair.second << "\n";
									homographymatrix_log << "relation: " << (relativePose_info.b_coplanar ? "homography" : "essential") << "\n";
									homographymatrix_log << "homography:"<< relativePose_homographyinfo.vec_inliers.size() << " in " << xI.cols() << "\n";
									for(int i = 0 ;i < 3; i++)
									{
										for(int j = 0 ;j < 3; j++)
										{
											homographymatrix_log << relativePose_homographyinfo.model_matrix(i,j)<<" ";
										}
										homographymatrix_log <<"\n";
									}
									homographymatrix_log << "relative pose:\n";
									Mat3 R_H = relativePose_homographyinfo.relativePose.rotation();
									Vec3 t_H = relativePose_homographyinfo.relativePose.translation();
									for (int i = 0; i < 3; i++)
									{
										for (int j = 0; j < 3; j++)
										{
											homographymatrix_log << R_H(i, j) << " ";
										}
										homographymatrix_log << t_H(i) << "\n";
									}
								}
								relative_poses.emplace(current_pair, pose_J);   //bc
								scoring_per_pair.emplace_back(scoring_angle, current_pair);
							}
							
						}
					}
				}
			}
		} // omp section
	}
	std::sort(scoring_per_pair.begin(), scoring_per_pair.end());
	// Since scoring is ordered in increasing order, reverse the order
	std::reverse(scoring_per_pair.begin(), scoring_per_pair.end());
	/////////bc debug start///////
	

	if (bdebug)
	{
		essentialmatrix_log.close();
		homographymatrix_log.close();
		std::ofstream initial_score_log;
		initial_score_log.open(stlplus::create_filespec(sOutDir, "initial_score_"+ file_suffix +".txt"));
		initial_score_log << "view i and view j :score\n";
		for (const auto& pair : scoring_per_pair)
		{
			initial_score_log << pair.second.first << " " << pair.second.second << ":   " << pair.first << "\n";
		}
		initial_score_log.close();
	}

	/////////bc debug end///////
	if (!scoring_per_pair.empty())
	{
		//initial_pair = scoring_per_pair.begin()->second;
		return true;
	}
	return false;

}

/// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
bool RobustMakeInitialPair3D(const Pair & current_pair,const SfM_Data& sfm_data_,const std::string sOut_directory_,
							const size_t max_iteration_count,const Matches_Provider* matches_provider_,Features_Provider* features_provider_,
							openMVG::tracks::SharedTrackVisibilityHelper* shared_track_visibility_helper_)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  const uint32_t
    I = std::min(current_pair.first, current_pair.second),
    J = std::max(current_pair.first, current_pair.second);

  if (sfm_data_.GetViews().count(I) == 0 ||
      sfm_data_.GetViews().count(J) == 0)
  {
    return false;
  }
  // a. Assert we have valid cameras
  const View
    * view_I = sfm_data_.GetViews().at(I).get(),
    * view_J = sfm_data_.GetViews().at(J).get();
  const Intrinsics::const_iterator
    iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
    iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

  if (iterIntrinsic_I == sfm_data_.GetIntrinsics().end() ||
      iterIntrinsic_J == sfm_data_.GetIntrinsics().end() )
  {
    return false;
  }

  const auto
    * cam_I = iterIntrinsic_I->second.get(),
    * cam_J = iterIntrinsic_J->second.get();
  if (!cam_I || !cam_J)
  {
    return false;
  }

  // b. Get common features between the two view
  // use the track to have a more dense match correspondence set
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);
///////////////////bc debug start/////////////////
bool bdebug = true;
if(bdebug)
{
  std::ofstream initial_tack_log(stlplus::create_filespec(sOut_directory_,"initial_track.txt"));
  initial_tack_log<<I<<"-"<<J<<"\n";
  std::set<Pair> set_initial_tacks;
  for(const auto& submaptrack_item:map_tracksCommon)
  {
    set_initial_tacks.insert(Pair(submaptrack_item.second.at(I),submaptrack_item.second.at(J)));
    //initial_tack_log<<submaptrack_item.second.at(I)<<" "<<submaptrack_item.second.at(J)<<"\n";
  }
  for(const auto& pair:set_initial_tacks)
  {
    initial_tack_log<<pair.first<<" "<<pair.second<<"\n";
  }
  initial_tack_log.close();
}

///////////////////bc debug end  /////////////////
  //-- Copy point to arrays
  const size_t n = map_tracksCommon.size();
  Mat xI(2,n), xJ(2,n);
  uint32_t cptIndex = 0;
  for (const auto & track_iter : map_tracksCommon)
  {
    auto iter = track_iter.second.cbegin();
    const uint32_t
      i = iter->second,
      j = (++iter)->second;

    Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
    xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
    feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
    xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
    ++cptIndex;
  }

  // c. Robust estimation of the relative pose
  RelativePose_MultiInfo relativePose_info;
  RelativePose_MultiInfo relativePose_essentialinfo;
  RelativePose_MultiInfo relativePose_homographyinfo;
  const std::pair<size_t, size_t>
    imageSize_I(cam_I->w(), cam_I->h()),
    imageSize_J(cam_J->w(), cam_J->h());
  if (bdebug)  //bc
  {
	  std::ofstream ransac_file(stlplus::create_filespec(sOut_directory_, "ransac_input.txt"));
	  ransac_file << "image size I(w,h):" << imageSize_I.first << " " << imageSize_I.second << "\n";
	  
	  ransac_file << "intrinsics I:\n";
	  Mat KI = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_I)->K();
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  ransac_file << KI(i, j) << " ";
		  }
		  ransac_file << "\n";
	  }
	  ransac_file << "image size J(w,h):" << imageSize_J.first << " " << imageSize_J.second << "\n";
	  ransac_file << "intrinsics J:\n";
	  Mat KJ = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_J)->K();
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  ransac_file << KJ(i, j) << " ";
		  }
		  ransac_file << "\n";
	  }
	  ransac_file << "xi yi xj yj:\n";
	  for (int i = 0; i < n; i++)
	  {
		  ransac_file << xI(0, i) << " " << xI(1, i) << " ";
		  ransac_file << xJ(0, i) << " " << xJ(1, i) << "\n";
	  }
	  ransac_file.close();
  }
  if (!robustRelativePose_MultiModel(
    cam_I, cam_J, xI, xJ, relativePose_info,imageSize_I, imageSize_J, max_iteration_count))
  {
    std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
      << std::endl;
    return false;
  }
  ////////bc debug start//////
  
  if (bdebug)  //bc
  {
	  std::ofstream initial_matrix_log;
	  initial_matrix_log.open(stlplus::create_filespec(sOut_directory_,"initial_matrix.txt"));
	  initial_matrix_log << "Pair " << view_I->id_view << " - " << view_J->id_view << "\n";
    initial_matrix_log <<"relation: "<<(relativePose_info.b_coplanar?"homography":"essential")<<"\n";
	  initial_matrix_log << "Essential:"<<relativePose_essentialinfo.vec_inliers.size()<<" in "<< xI .cols()<<"\n";
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  initial_matrix_log << relativePose_essentialinfo.model_matrix(i,j)<<" ";
		  }
		  initial_matrix_log << "\n";
	  }
	  initial_matrix_log << "relativePose:\n";
	  Mat3 R_E = relativePose_essentialinfo.relativePose.rotation();
		Vec3 t_E = relativePose_essentialinfo.relativePose.translation();
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  initial_matrix_log << R_E(i, j) << " ";
		  }
		  initial_matrix_log << t_E(i) << "\n";
	  }


    initial_matrix_log << "homography:"<< relativePose_homographyinfo.vec_inliers.size() << " in " << xI.cols() << "\n";
    for(int i = 0 ;i < 3; i++)
    {
        for(int j = 0 ;j < 3; j++)
        {
          initial_matrix_log << relativePose_homographyinfo.model_matrix(i,j)<<" ";
        }
        initial_matrix_log <<"\n";
    }
    initial_matrix_log << "relative pose:\n";
    Mat3 R_H = relativePose_homographyinfo.relativePose.rotation();
    Vec3 t_H = relativePose_homographyinfo.relativePose.translation();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
          initial_matrix_log << R_H(i, j) << " ";
        }
        initial_matrix_log << t_H(i) << "\n";
    }


	  initial_matrix_log << "inliers:\n";
	  for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
	  {
		  openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
		  std::advance(iterT, inlier_idx);
		  tracks::submapTrack::const_iterator iter = iterT->second.begin();
		  assert(map_tracksCommon.at(inlier_idx).at(I) == iter->second);
		  assert(map_tracksCommon.at(inlier_idx).at(J) == (++iter)->second);
		  const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
		  const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
		  initial_matrix_log << iter->second << " " << featI(0) << " " << featI(1) << " ";
		  initial_matrix_log << (++iter)->second << " " << featJ(0) << " " << featJ(1) << "\n";
	  }
	  initial_matrix_log.close();
  }
  ////////bc debug end//////
  std::cout << "A-Contrario initial pair residual: "
    << relativePose_info.found_residual_precision << std::endl;
//   // Bound min precision at 1 pix.
//   relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

//   const bool bRefine_using_BA = true;
//   if (bRefine_using_BA)
//   {
//     // Refine the defined scene
//     SfM_Data tiny_scene;
//     tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
//     tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
//     tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
//     tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

//     // Init poses
//     const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
//     const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

//     // Init structure
//     Landmarks & landmarks = tiny_scene.structure;

//     for (const auto & track_iterator : map_tracksCommon)
//     {
//       // Get corresponding points
//       auto iter = track_iterator.second.cbegin();
//       const uint32_t
//         i = iter->second,
//         j = (++iter)->second;

//       const Vec2
//         x1 = features_provider_->feats_per_view[I][i].coords().cast<double>(),
//         x2 = features_provider_->feats_per_view[J][j].coords().cast<double>();

//       Vec3 X;
//       if (Triangulate2View(
//             Pose_I.rotation(),
//             Pose_I.translation(),
//             (*cam_I)(cam_I->get_ud_pixel(x1)),
//             Pose_J.rotation(),
//             Pose_J.translation(),
//             (*cam_J)(cam_J->get_ud_pixel(x2)),
//             X,
//             triangulation_method_))
//       {
//         Observations obs;
//         obs[view_I->id_view] = Observation(x1, i);
//         obs[view_J->id_view] = Observation(x2, j);
//         landmarks[track_iterator.first].obs = std::move(obs);
//         landmarks[track_iterator.first].X = X;
//       }
//     }
//     Save(tiny_scene, stlplus::create_filespec(sOut_directory_, "initialPair.ply"), ESfM_Data(ALL));

//     // - refine only Structure and Rotations & translations (keep intrinsic constant)
//     Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
//     options.linear_solver_type_ = ceres::DENSE_SCHUR;
//     Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
//     if (!bundle_adjustment_obj.Adjust(tiny_scene,
//         Optimize_Options
//         (
//           Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
//           Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
//           Structure_Parameter_Type::ADJUST_ALL) // Adjust structure
//         )
//       )
//     {
//       return false;
//     }

//     // Save computed data
//     const Pose3 pose_I = sfm_data_.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
//     const Pose3 pose_J = sfm_data_.poses[view_J->id_pose] = tiny_scene.poses[view_J->id_pose];
//     map_ACThreshold_.insert({I, relativePose_info.found_residual_precision});
//     map_ACThreshold_.insert({J, relativePose_info.found_residual_precision});
//     set_remaining_view_id_.erase(view_I->id_view);
//     set_remaining_view_id_.erase(view_J->id_view);

//     // List inliers and save them
//     std::ofstream initial_triangulation_log;////bc
//     if(bdebug)
//     {
// 		std::ofstream ba_initial_essential_log(stlplus::create_filespec(sOut_directory_, "ba_initial_essentials.txt"));
// 		ba_initial_essential_log << "ba relative pose:\n";
// 		Pose3 relative_poseji = pose_J * pose_I.inverse();
// 		Mat3 relative_rji = relative_poseji.rotation();
// 		Vec3 relative_tji = relative_poseji.translation();
// 		for (int i = 0; i < 3; i++)
// 		{
// 			for (int j = 0; j < 3; j++)
// 			{
// 				ba_initial_essential_log << relative_rji(i, j) << " ";
// 			}
// 			ba_initial_essential_log << relative_tji(i) << "\n";
// 		}
// 		Mat3 tji_ni;
// 		tji_ni << 0,                -relative_tji(2), relative_tji(1),
// 			      relative_tji(2),  0,                -relative_tji(0),
// 			      -relative_tji(1), relative_tji(0),  0;
// 		Mat3 essential_matrix = tji_ni * relative_rji;
// 		ba_initial_essential_log << "ba essential matrix:\n";
// 		for (int i = 0; i < 3; i++)
// 		{
// 			for (int j = 0; j < 3; j++)
// 			{
// 				ba_initial_essential_log << essential_matrix(i, j) << " ";
// 			}
// 			ba_initial_essential_log << "\n";
// 		}
// 		ba_initial_essential_log << "essential_error\n";
// 		std::vector<double> essential_error;
// 		for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
// 		{
// 			openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
// 			std::advance(iterT, inlier_idx);
// 			tracks::submapTrack::const_iterator iter = iterT->second.begin();
// 			assert(map_tracksCommon.at(inlier_idx).at(I) == iter->second);
// 			assert(map_tracksCommon.at(inlier_idx).at(J) == (++iter)->second);
// 			const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
// 			const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
// 			Vec2 norm_featI=cam_I->ima2cam(featI);
// 			Vec2 norm_featJ=cam_J->ima2cam(featJ);
// 			Vec3 norm_featI_3, norm_featJ_3;
// 			norm_featI_3 << norm_featI(0), norm_featI(1), 1.0;
// 			norm_featJ_3 << norm_featJ(0), norm_featJ(1), 1.0;
// 			Vec error = norm_featJ_3.transpose() * essential_matrix * norm_featI_3;
// 			essential_error.push_back(error(0));
			
// 		}
// 		std::sort(essential_error.begin(), essential_error.end());
// 		for (size_t i = 0; i < essential_error.size(); i++)
// 		{
// 			ba_initial_essential_log << essential_error[i] << "\n";
// 		}
// 		ba_initial_essential_log.close();
//       initial_triangulation_log.open(stlplus::create_filespec(sOut_directory_,"initial_triangulation_status.txt"));
// 	  initial_triangulation_log << view_I->id_view << "-" << view_J->id_view << "\n";
// 	}
//     for (const auto & landmark_entry : tiny_scene.GetLandmarks())
//     {
//       const IndexT trackId = landmark_entry.first;
//       const Landmark & landmark = landmark_entry.second;
//       const Observations & obs = landmark.obs;
//       Observations::const_iterator
//         iterObs_xI = obs.find(view_I->id_view),
//         iterObs_xJ = obs.find(view_J->id_view);

//       const Observation & ob_xI = iterObs_xI->second;
//       const Observation & ob_xJ = iterObs_xJ->second;
//       const Vec2
//         ob_xI_ud = cam_I->get_ud_pixel(ob_xI.x),
//         ob_xJ_ud = cam_J->get_ud_pixel(ob_xJ.x);

//       const double angle = AngleBetweenRay(
//         pose_I, cam_I, pose_J, cam_J, ob_xI_ud, ob_xJ_ud);
//       const Vec2 residual_I = cam_I->residual(pose_I(landmark.X), ob_xI.x);
//       const Vec2 residual_J = cam_J->residual(pose_J(landmark.X), ob_xJ.x);
//       if (angle > 2.0 &&
//           CheiralityTest((*cam_I)(ob_xI_ud), pose_I,
//                          (*cam_J)(ob_xJ_ud), pose_J,
//                          landmark.X) &&
//           residual_I.norm() < relativePose_info.found_residual_precision &&
//           residual_J.norm() < relativePose_info.found_residual_precision)
//       {
//         sfm_data_.structure[trackId] = landmarks[trackId];
// 		if (bdebug)
// 		{
// 			initial_triangulation_log << ob_xI.id_feat << " " << ob_xJ.id_feat 
// 				<<"("<< ob_xI.x(0)<<" "<< ob_xI.x(1)<<" "<< ob_xJ.x(0)<<" "<< ob_xJ.x(1)<<")"
// 				<< ": triangulated\n";  //bc
// 		}
//       }
//       else if(bdebug)  //bc
//       {
//         initial_triangulation_log << ob_xI.id_feat<<" "<<ob_xJ.id_feat<<": ";  //bc
//         if(!(angle > 2.0))
//         {
//           initial_triangulation_log<<" angle <= 2.0 ";
//         }
//         else if(!CheiralityTest((*cam_I)(ob_xI_ud), pose_I,
//                          (*cam_J)(ob_xJ_ud), pose_J,
//                          landmark.X))
//         {
//           initial_triangulation_log<<" CheiralityTest failed ";
//         }
//         else if(!(residual_I.norm() < relativePose_info.found_residual_precision ))
//         {
//           initial_triangulation_log<<" Reprjection of "<<ob_xI.id_feat<<" is too large ";
//         }
//         else if(!(residual_J.norm() < relativePose_info.found_residual_precision))
//         {
//           initial_triangulation_log<<" Reprjection of "<<ob_xJ.id_feat<<" is too large ";
//         }
//         initial_triangulation_log<<"\n";
//       }
      
//     }
//     initial_triangulation_log.close();
//     // Save outlier residual information
//     Histogram<double> histoResiduals;
//     std::cout << "\n"
//       << "=========================\n"
//       << " MSE Residual InitialPair Inlier:\n";
//     ComputeResidualsHistogram(&histoResiduals);
//     std::cout << "=========================" << std::endl;

//     if (!sLogging_file_.empty())
//     {
//       using namespace htmlDocument;
//       html_doc_stream_->pushInfo(htmlMarkup("h1","Essential Matrix."));
//       std::ostringstream os;
//       os << std::endl
//         << "-------------------------------" << "<br>"
//         << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
//         << view_I->s_Img_path << ","
//         << view_J->s_Img_path << "<br>"
//         << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
//         << "-- Resection status: " << "OK" << "<br>"
//         << "-- Nb points used for robust Essential matrix estimation: "
//         << xI.cols() << "<br>"
//         << "-- Nb points validated by robust estimation: "
//         << sfm_data_.structure.size() << "<br>"
//         << "-- % points validated: "
//         << sfm_data_.structure.size()/static_cast<float>(xI.cols())
//         << "<br>"
//         << "-------------------------------" << "<br>";
//       html_doc_stream_->pushInfo(os.str());

//       html_doc_stream_->pushInfo(htmlMarkup("h2",
//         "Residual of the robust estimation (Initial triangulation). Thresholded at: "
//         + toString(relativePose_info.found_residual_precision)));

//       html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of residuals"));

//       const std::vector<double> xBin = histoResiduals.GetXbinsValue();
//       const auto range = autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

//       htmlDocument::JSXGraphWrapper jsxGraph;
//       jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
//       jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
//       jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
//         relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
//       jsxGraph.UnsuspendUpdate();
//       jsxGraph.setViewport(range);
//       jsxGraph.close();
//       html_doc_stream_->pushInfo(jsxGraph.toStr());

//       html_doc_stream_->pushInfo("<hr>");

//       std::ofstream htmlFileStream( std::string(stlplus::folder_append_separator(sOut_directory_) +
//         "Reconstruction_Report.html").c_str());
//       htmlFileStream << html_doc_stream_->getDoc();
//     }
//   }
  return !sfm_data_.structure.empty();
}

bool Choose_InitialPair(const SfM_Data& sfm_data, const std::set<IndexT>& registered_views,
	const Matches_Provider* matches_provider, const Matches_Provider* notablematches_provider, Features_Provider* features_provider,
	const Pair_Set& tried_initialpairs,
	const std::string sOutDir, const int& iteration_i, std::vector<Pair>& initial_pairs,const size_t max_iteration_count)
{
	//InitLandmarkTracks
	openMVG::tracks::STLMAPTracks map_tracks;
	std::shared_ptr<openMVG::tracks::SharedTrackVisibilityHelper> shared_track_visibility_helper = InitLandmarkTracks(matches_provider, map_tracks);

	const double MaxAngleError = 2;
	const double MaxresidualError = 0.08;
	openMVG::system::Timer timer;
	Hash_Map<Pair, Pose3> relative_poses;
	std::vector<std::pair<double, Pair>> scoring_per_pair;
	if (initial_pairs.empty())
	{
	FindInitialImagePairs(sfm_data, registered_views,
		matches_provider, features_provider, tried_initialpairs,
		relative_poses, scoring_per_pair,
		sOutDir, iteration_i,"",max_iteration_count,shared_track_visibility_helper.get());
		std::cout << std::endl << "find_initialpair took (s): " << timer.elapsed() << std::endl;
		std::cout << "choose initial Pair " << scoring_per_pair.begin()->second.first << " - " << scoring_per_pair.begin()->second.second << "\n";
	
		initial_pairs.push_back(scoring_per_pair.begin()->second);
	}
	RobustMakeInitialPair3D(*initial_pairs.begin(),sfm_data,sOutDir,
							max_iteration_count,matches_provider, features_provider,shared_track_visibility_helper.get());
	return true;

	Hash_Map<Pair, Pose3> notablerelative_poses;
	std::vector<std::pair<double, Pair>> notablescoring_per_pair;
	FindInitialImagePairs(sfm_data, registered_views,
		notablematches_provider, features_provider, tried_initialpairs,
		notablerelative_poses, notablescoring_per_pair,
		sOutDir, iteration_i,"notable", max_iteration_count, shared_track_visibility_helper.get());
	std::vector<double> vec_residualErrors;
	std::vector<double> vec_angularErrors;
	bool bdebug = true;

	std::ofstream errorfile;
	if (bdebug)
	{
		errorfile.open(stlplus::create_filespec(sOutDir, "errorfile.txt"));
		errorfile << "view1 view2 angularErrorDegree dCenterResidual score score_notable\n";
	}

	//bc debug start//
	std::map<Pair, double> map_notablescore;
	for (const auto& pair_item : notablescoring_per_pair)
	{
		map_notablescore.emplace(pair_item.second, pair_item.first);
	}
	//bc debug end//
	for (std::vector<std::pair<double, Pair>>::iterator iter = scoring_per_pair.begin();
		iter != scoring_per_pair.end() ; iter++)
	{
		Pair pair = iter->second;
		if (!relative_poses.count(pair) || !notablerelative_poses.count(pair))
		{
			continue;
		}
		Pose3 relative_pose = relative_poses.at(pair);
		
		Pose3 notablerelative_pose = notablerelative_poses.at(pair);
		// Compute statistics and export them
		// -a. distance between camera center
		// -b. angle between rotation matrix

		// -a. distance between camera center
		
		
		const double dResidual = (relative_pose.center() - notablerelative_pose.center()).norm();
		vec_residualErrors.push_back(dResidual);

		

		// -b. angle between rotation matrix
		
		
		const Mat3 R1 = relative_pose.rotation(); 
		const Mat3 R2 = notablerelative_pose.rotation(); 

		const double angularErrorDegree = R2D(getRotationMagnitude(R1 * R2.transpose()));
		vec_angularErrors.push_back(angularErrorDegree);
		if (dResidual < MaxresidualError&& angularErrorDegree < MaxAngleError)
		{
			initial_pairs.push_back(pair);
		}
		if (bdebug)
		{
			errorfile << pair.first << " " << pair.second << ": " << angularErrorDegree << " " << dResidual
				<< " " << iter->first << " " << map_notablescore.at(iter->second)
				<< "\n";
		}

	}
	//bc debug start//
	if (bdebug)
	{
		errorfile.close();
		std::ofstream rot_errorfile(stlplus::create_filespec(sOutDir, "rot_errorfile.txt"));
		std::sort(vec_angularErrors.begin(), vec_angularErrors.end());
		utils::AllMinMaxMedianMean(rot_errorfile, vec_angularErrors);
		rot_errorfile.close();

		std::ofstream center_errorfile(stlplus::create_filespec(sOutDir, "center_errorfile.txt"));
		std::sort(vec_residualErrors.begin(), vec_residualErrors.end());
		utils::AllMinMaxMedianMean(center_errorfile, vec_residualErrors);
		center_errorfile.close();
	}
	

	//bc debug end//
	return true;


}



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
  std::string notable_features_dir = "";
  std::pair<std::string,std::string> initialPairString("","");
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
  bool b_use_motion_priors = false;
  int triangulation_method = static_cast<int>(ETriangulationMethod::DEFAULT);
  size_t initial_max_iteration_count = 256;

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
  cmd.add( make_option('n', notable_features_dir, "notable_features_dir"));
  cmd.add( make_option('d', initial_max_iteration_count, "initial_max_iteration_count"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n[Optional]\n"
		<< "[-n|--notable_features_dir] directory where the notable features stored\n"
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
  std::cout<<"initial_max_iteration_count: "<<initial_max_iteration_count <<"\n";
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

  // notable features reading
  std::cout << "notable features reading\n";
  bool bContinue = true;
  std::map<IndexT, std::set<IndexT>> notable_features;
  for (Views::const_iterator iter = raw_sfm_data.GetViews().begin();
	  iter != raw_sfm_data.GetViews().end() && bContinue; ++iter)
  {
	  const IndexT  view_id = iter->first;
	  {
		  const std::string sImageName = stlplus::create_filespec(raw_sfm_data.s_root_path, iter->second->s_Img_path);
		  const std::string basename = stlplus::basename_part(sImageName);
		  const std::string featFile = stlplus::create_filespec(notable_features_dir, basename, ".txt");

		  
		  if (!stlplus::file_exists(featFile) )
		  {
			  //std::cerr << "Invalid feature files for the view: " << sImageName << std::endl;
			  //bContinue = false;
			  continue;
		  }
		  
		  const features::PointFeatures & pointfeatures = feats_provider->getFeatures(view_id);
		  std::ifstream file(featFile);
		  IndexT feat_id;
		  float x, y;
		  size_t num_features = 0;
		  while (file >> feat_id >> x >> y)
		  {
			  if (pointfeatures[feat_id].x() != x || pointfeatures[feat_id].y() != y)
			  {
				  std::cerr << "Error:the feature is inconsistent.\n";
				  bContinue = false;
			  }
			  else
			  {
				  num_features++;
				  notable_features[view_id].insert(feat_id);
			  }
			  if (!bContinue) break;
			  
		  }
		  std::cout << "loaded " << num_features << " for view"<< view_id <<" from " << basename << ".txt\n";
		  
	  }
  }
  if (!bContinue)
  {
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
  //notable matches filtering
  std::cout << "notable matches filtering\n";
  std::shared_ptr<Matches_Provider> notablematches_provider = std::make_shared<Matches_Provider>();
  for (const auto& pairwisematch : matches_provider->pairWise_matches_)
  {
	  const IndexT I = pairwisematch.first.first;
	  const IndexT J = pairwisematch.first.second;
	  if (!notable_features.count(I) || !notable_features.count(J))  //for last 2 frame ,there are none notable features
	  {
		  continue;
	  }
	  matching::IndMatches notable_matches;
	  for (const auto& indmatch : pairwisematch.second)
	  {
		  if (notable_features.at(I).count(indmatch.i_) && notable_features.at(J).count(indmatch.j_))
		  {
			  notable_matches.emplace_back(indmatch);
		  }
	  }
	  if (!notable_matches.empty())
	  {
		  notablematches_provider->pairWise_matches_.insert({ { I,J }, std::move(notable_matches) });
	  }
  }
  // check output directory
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

  openMVG::system::Timer timer;
  std::set<IndexT> registered_views;
  Pair_Set initial_tried_pairs;
  std::set<IndexT> first_initial_tried_images;
  Pair_Set tried_initialpairs;
  int iteration_i = 1;
  size_t pose_before = 0;
  const size_t num_views = raw_sfm_data.GetViews().size();
  size_t num_reconstructions = 0;
  while(num_views - registered_views.size() > 3)
  {
    std::cout<<"/////////////////////////\n"
             <<"///////iteration "<<iteration_i<<"///////\n"
             <<"/////////////////////////\n";
    
    // find the initial pairs in not yet registered views
    std::vector<Pair> initial_pairs;
    
    {
	if (!initialPairString.first.empty() && !initialPairString.second.empty())
	{
		Pair initialPairIndex;
		if (!computeIndexFromImageNames(raw_sfm_data, initialPairString, initialPairIndex))
		{
			std::cerr << "Could not find the initial pairs <" << initialPairString.first
			<<  ", " << initialPairString.second << ">!\n";
		return EXIT_FAILURE;
		}
		initial_pairs.push_back(initialPairIndex);
	}
      if(!Choose_InitialPair(raw_sfm_data,registered_views,
                                matches_provider.get(),notablematches_provider.get(),
								feats_provider.get(),tried_initialpairs, 
                                sOutDir, iteration_i,initial_pairs,initial_max_iteration_count))
      {
        std::cout<<"There are not good initial pairs for next reconstruction\n";
        break;
      }
	  
    }
	break;   ///bc

	std::string iteration_sOutDir = stlplus::create_filespec(sOutDir, std::to_string(iteration_i));
	if (!stlplus::folder_exists(iteration_sOutDir))
	{
		if (!stlplus::folder_create(iteration_sOutDir))
		{
			std::cerr << "\nCannot create the output directory" << std::endl;
		}
	}

    //try to reconstruction with initial pair
    //  * If reconstruct successfully,exit loop and find next initial pairs for not yet registered views
	bool brecons = false;
    for(const Pair& initial_pair:initial_pairs)
    {
      openMVG::system::Timer re_timer;
      SequentialSfMReconstructionEngine sfmEngine(
        raw_sfm_data,
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
	  bool bProcess = sfmEngine.Process();
      if (bProcess && sfmEngine.Get_SfM_Data().GetPoses().size() > 3)
      {
		  brecons = true;
        num_reconstructions ++;
        std::cout << "...Generating SfM_Report.html" << std::endl;
        Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
          stlplus::create_filespec(iteration_sOutDir, "SfMReconstruction_Report.html"));

        //-- Export to disk computed scene (data & visualizable results)
        std::cout << "...Export SfM_Data to disk." << std::endl;
        Save(sfmEngine.Get_SfM_Data(),
          stlplus::create_filespec(iteration_sOutDir, "sfm_data", ".bin"),
          ESfM_Data(ALL));

        Save(sfmEngine.Get_SfM_Data(),
          stlplus::create_filespec(iteration_sOutDir, "cloud_and_poses", ".ply"),
          ESfM_Data(ALL));
        if(raw_sfm_data.GetPoses().size() > 0)
        {
          std::cout<<"Error:raw_sfm_data has been modified\n";
          return EXIT_SUCCESS;
        }

        //update the views registered into `registered_views`
        const SfM_Data& iteration_sfm_data = sfmEngine.Get_SfM_Data();
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
	break;
  }
  std::cout<<"\n////////////Holistic Static//////////////\n";
  std::cout<<"#Images registered:"<<registered_views.size()<<"\n";
  std::cout<<"#Images not yet registered:"<<num_views - registered_views.size()<<"\n";
  std::cout << std::endl << " Total Ac-Sfm took (s): " <<timer.elapsed() << std::endl;
  std::cout<<"/////////////////////\n";

  return EXIT_SUCCESS;
}

