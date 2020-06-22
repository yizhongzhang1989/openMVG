// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG_IMU/sfm/pipelines/sequential/sequential_SfM_robust_initialization.hpp"
#include "openMVG_IMU/sfm/pipelines/sfm_robust_multimodel_estimation.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#include <ceres/types.h>
#include <functional>
#include <iostream>
#include <utility>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;

SequentialSfMReconstructionEngine_Robust_Initialization::SequentialSfMReconstructionEngine_Robust_Initialization(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile,
  const size_t initial_max_iteration_count)
  : SequentialSfMReconstructionEngine_General(sfm_data, soutDirectory,sloggingFile),
  initial_max_iteration_count_(initial_max_iteration_count)
{
  b_robust_initialization_of_imu_ = false;
}

SequentialSfMReconstructionEngine_Robust_Initialization::~SequentialSfMReconstructionEngine_Robust_Initialization()
{
  
}


bool SequentialSfMReconstructionEngine_Robust_Initialization::Process_Robust_Initialization() {

  //-------------------
  //-- Incremental reconstruction
  //-------------------
  /////////BC DEBUG START////////
	bool bdebug = true;  //bc
	std::ofstream intrinsic_log;
	if (bdebug)
	{
		intrinsic_log.open(stlplus::create_filespec(sOut_directory_, "intrinsic_log.txt"));
		intrinsic_log << "input intrinsic: \n";
		Mat3 K = dynamic_cast<const cameras::Pinhole_Intrinsic*>(sfm_data_.GetIntrinsics().at(0).get())->K();
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				intrinsic_log << K(i, j) << " ";
			}
			intrinsic_log << "\n";
		}
	}
  /////////BC DEBUG END////////
  if (!InitLandmarkTracks())
    return false;
  
  // Initial pair choice
  openMVG::system::Timer initial_timer;
  if (initial_pair_ == Pair(0,0))
  {
    if (!RobustAutomaticInitialPairChoice(initial_pair_))
    {
      intrinsic_log.close();
      bool flag = Process_KnownTracks();   //turn the level of initialization down 
      return flag;
    }
  }
  std::cout << std::endl << " Total Initialization took (s): " << initial_timer.elapsed() << std::endl;
  // Else a starting pair was already initialized before

  // Initial pair Essential Matrix and [R|t] estimation.
  if (!RobustMakeInitialPair3D(initial_pair_))
    return false;
  //return true;  //must deleted
  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  size_t resectionGroupIndex = 0;
  std::vector<uint32_t> vec_possible_resection_indexes;
  
  std::vector<uint32_t> registration_order;  //bc
  if (bdebug)
  {
	  intrinsic_log << "initial intrinsic"<<"("<< initial_pair_ .first<<","<< initial_pair_ .second<<"): \n";
	  Mat3 K = dynamic_cast<const cameras::Pinhole_Intrinsic*>(sfm_data_.GetIntrinsics().at(0).get())->K();
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  intrinsic_log << K(i, j) << " ";
		  }
		  intrinsic_log << "\n";
	  }
	  Save(sfm_data_,
		  stlplus::create_filespec(sOut_directory_, "sfm_data_" + std::to_string(registration_order.size()) + "_" 
									+ std::to_string(initial_pair_.first)+std::to_string(initial_pair_.second), ".bin"),
		  ESfM_Data(ALL));
  }
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (const auto & iter : vec_possible_resection_indexes)
    {
      bImageAdded |= Resection(iter);
      set_remaining_view_id_.erase(iter);
	  if (bdebug) //bc
	  {
		  registration_order.push_back(iter);
		  Save(sfm_data_,
			  stlplus::create_filespec(sOut_directory_, "sfm_data_"+std::to_string(registration_order.size())+"_"+std::to_string(iter), ".bin"),
			  ESfM_Data(ALL));
		               
	  }
    }

    if (bImageAdded)
    {
      // Scene logging as ply for visual debug
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
      Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

      // Perform BA until all point are under the given precision
      do
      {
        BundleAdjustment();
      }
      while (badTrackRejector(4.0, 50));
      eraseUnstablePosesAndObservations(sfm_data_);
    }
	if (bdebug) //bc
	{
		intrinsic_log <<" Iteration "<< registration_order.size() <<" intrinsic: \n";
		Mat3 K = dynamic_cast<const cameras::Pinhole_Intrinsic*>(sfm_data_.GetIntrinsics().at(0).get())->K();
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				intrinsic_log << K(i, j) << " ";
			}
			intrinsic_log << "\n";
		}
		Save(sfm_data_,
			stlplus::create_filespec(sOut_directory_, "ba_sfm_data_" + std::to_string(registration_order.size()), ".bin"),
			ESfM_Data(ALL));
	}
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  if (badTrackRejector(4.0, 0))
  {
    eraseUnstablePosesAndObservations(sfm_data_);
  }
  /////////BC START  DEBUG////////
  if (bdebug)
  {
	  intrinsic_log.close();
	  std::ofstream registration_order_log(stlplus::create_filespec(sOut_directory_, "registration_order.txt"));
	  for (const auto& view_id : registration_order)
	  {
		  registration_order_log << view_id << " ";
	  }
	  registration_order_log.close();
  }
  
  /////////BC END  DEBUG////////
  //-- Reconstruction done.
  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
    << " from " << sfm_data_.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  std::cout << "\nHistogram of residuals:\n" << h.ToString() << std::endl;

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion process finished.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Structure from Motion (statistics):<br>"
      << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
      << " from " <<sfm_data_.GetViews().size() << " input images.<br>"
      << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());

    html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

    const std::vector<double> xBin = h.GetXbinsValue();
    const auto range = autoJSXGraphViewport<double>(xBin, h.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("3DtoImageResiduals",600,300);
    jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());
  }
  return true;
}

bool SequentialSfMReconstructionEngine_Robust_Initialization::Process_KnownTracks() {

  //-------------------
  //-- Incremental reconstruction
  //-------------------
  /////////BC DEBUG START////////
	bool bdebug = true;  //bc
  std::cout<<"Initialization of OpenMVG\n";
  std::ofstream intrinsic_log;
	if (bdebug)
	{
		intrinsic_log.open(stlplus::create_filespec(sOut_directory_, "intrinsic_log.txt"));
		intrinsic_log << "input intrinsic: \n";
		Mat3 K = dynamic_cast<const cameras::Pinhole_Intrinsic*>(sfm_data_.GetIntrinsics().at(0).get())->K();
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				intrinsic_log << K(i, j) << " ";
			}
			intrinsic_log << "\n";
		}
	}
  /////////BC DEBUG END////////
  // Initial pair choice
  if (initial_pair_ == Pair(0,0))
  {
    if (!AutomaticInitialPairChoice(initial_pair_))
    {
      // Cannot find a valid initial pair, try to set it by hand?
      if (!ChooseInitialPair(initial_pair_))
      {
        return false;
      }
    }
  }
  // Else a starting pair was already initialized before

  // Initial pair Essential Matrix and [R|t] estimation.
  if (!MakeInitialPair3D(initial_pair_))
    return false;
  //return true;  //must deleted
  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  size_t resectionGroupIndex = 0;
  std::vector<uint32_t> vec_possible_resection_indexes;
  
  std::vector<uint32_t> registration_order;  //bc
  if (bdebug)
  {
    intrinsic_log << "initial intrinsic"<<"("<< initial_pair_ .first<<","<< initial_pair_ .second<<"): \n";
	  Mat3 K = dynamic_cast<const cameras::Pinhole_Intrinsic*>(sfm_data_.GetIntrinsics().at(0).get())->K();
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  intrinsic_log << K(i, j) << " ";
		  }
		  intrinsic_log << "\n";
	  }
	  Save(sfm_data_,
		  stlplus::create_filespec(sOut_directory_, "sfm_data_" + std::to_string(registration_order.size()) + "_" 
									+ std::to_string(initial_pair_.first)+std::to_string(initial_pair_.second), ".bin"),
		  ESfM_Data(ALL));
  }
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (const auto & iter : vec_possible_resection_indexes)
    {
      bImageAdded |= Resection(iter);
      set_remaining_view_id_.erase(iter);
	  if (bdebug) //bc
	  {
		  registration_order.push_back(iter);
		  Save(sfm_data_,
			  stlplus::create_filespec(sOut_directory_, "sfm_data_"+std::to_string(registration_order.size())+"_"+std::to_string(iter), ".bin"),
			  ESfM_Data(ALL));
		               
	  }
    }

    if (bImageAdded)
    {
      // Scene logging as ply for visual debug
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
      Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

      // Perform BA until all point are under the given precision
      do
      {
        BundleAdjustment();
      }
      while (badTrackRejector(4.0, 50));
      eraseUnstablePosesAndObservations(sfm_data_);
    }
	if (bdebug) //bc
	{
		intrinsic_log <<" Iteration "<< registration_order.size() <<" intrinsic: \n";
		Mat3 K = dynamic_cast<const cameras::Pinhole_Intrinsic*>(sfm_data_.GetIntrinsics().at(0).get())->K();
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				intrinsic_log << K(i, j) << " ";
			}
			intrinsic_log << "\n";
		}
		Save(sfm_data_,
			stlplus::create_filespec(sOut_directory_, "ba_sfm_data_" + std::to_string(registration_order.size()), ".bin"),
			ESfM_Data(ALL));
	}
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  if (badTrackRejector(4.0, 0))
  {
    eraseUnstablePosesAndObservations(sfm_data_);
  }
  /////////BC START  DEBUG////////
  if (bdebug)
  {
	  intrinsic_log.close();
	  std::ofstream registration_order_log(stlplus::create_filespec(sOut_directory_, "registration_order.txt"));
	  for (const auto& view_id : registration_order)
	  {
		  registration_order_log << view_id << " ";
	  }
	  registration_order_log.close();
  }
  
  /////////BC END  DEBUG////////
  //-- Reconstruction done.
  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
    << " from " << sfm_data_.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  std::cout << "\nHistogram of residuals:\n" << h.ToString() << std::endl;

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion process finished.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Structure from Motion (statistics):<br>"
      << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
      << " from " <<sfm_data_.GetViews().size() << " input images.<br>"
      << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());

    html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

    const std::vector<double> xBin = h.GetXbinsValue();
    const auto range = autoJSXGraphViewport<double>(xBin, h.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("3DtoImageResiduals",600,300);
    jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());
  }
  return true;
}



bool SequentialSfMReconstructionEngine_Robust_Initialization::RobustAutomaticInitialPairChoice(Pair & initial_pair) const
{
  // select a pair that have the largest baseline (mean angle between its bearing vectors).

  const unsigned iMin_inliers_count = 100;
  const float fRequired_min_angle = 3.0f;
  const float fLimit_max_angle = 60.0f; // More than 60 degree, we cannot rely on matches for initial pair seeding

  // List Views that support valid intrinsic (view that could be used for Essential matrix computation)
  std::set<IndexT> valid_views;
  for (Views::const_iterator it = sfm_data_.GetViews().begin();
    it != sfm_data_.GetViews().end(); ++it)
  {
    const View * v = it->second.get();
    if (sfm_data_.GetIntrinsics().count(v->id_intrinsic))
      valid_views.insert(v->id_view);
  }

  if (valid_views.size() < 2)
  {
    return false; // There is not view that support valid intrinsic data
  }

  std::vector<std::pair<double, Pair>> scoring_per_pair;
  //////////BC DEBUG START///////
	bool bdebug = true;
	std::ofstream essentialmatrix_log;
	std::ofstream homographymatrix_log;
	if(bdebug)
	{
		essentialmatrix_log.open(stlplus::create_filespec(sOut_directory_,"total_initial_essential_matrix.txt"));
		homographymatrix_log.open(stlplus::create_filespec(sOut_directory_, "total_initial_homography_matrix.txt"));
	}
	//////////BC DEBUG END///////
  // Compute the relative pose & the 'baseline score'
  C_Progress_display my_progress_bar( matches_provider_->pairWise_matches_.size(),
    std::cout,
    "Robust automatic selection of an initial pair:\n" );
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (const std::pair<Pair, IndMatches> & match_pair : matches_provider_->pairWise_matches_)
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
          * view_I = sfm_data_.GetViews().at(I).get(),
          * view_J = sfm_data_.GetViews().at(J).get();
        const Intrinsics::const_iterator
          iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
          iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

        const auto
          cam_I = iterIntrinsic_I->second.get(),
          cam_J = iterIntrinsic_J->second.get();
        if (cam_I && cam_J)
        {
          openMVG::tracks::STLMAPTracks map_tracksCommon;
          shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

          // Copy points correspondences to arrays for relative pose estimation
          const size_t n = map_tracksCommon.size();
          Mat xI(2,n), xJ(2,n);
          size_t cptIndex = 0;
          for (const auto & track_iter : map_tracksCommon)
          {
            auto iter = track_iter.second.cbegin();
            const uint32_t i = iter->second;
            const uint32_t j = (++iter)->second;

            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
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
          bool initialization_flag = false;
          if (!b_robust_initialization_of_imu_)
          {
          if (robustRelativePose_MultiModel(
                cam_I, cam_J,
                xI, xJ, relativePose_info,
                relativePose_essentialinfo,
						    relativePose_homographyinfo,
                {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
                initial_max_iteration_count_)
              && relativePose_info.vec_inliers.size() > iMin_inliers_count)
              {
                initialization_flag = true;
              }
          }
          else
          {
              const Pose3 imu_pose_I = imu_data_.poses.at(I);
              const Pose3 imu_pose_J = imu_data_.poses.at(J);
              Pose3 imu_relative_pose = imu_pose_J * imu_pose_I.inverse();
              if(robustRelativePose_IMU(
                cam_I, cam_J,
                xI, xJ, relativePose_info,imu_relative_pose,
                {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
                initial_max_iteration_count_)
              && relativePose_info.vec_inliers.size() > iMin_inliers_count)
              {
                  relativePose_essentialinfo = relativePose_info;
                  initialization_flag = true;
              }
          }
              
          if(initialization_flag)
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
              const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
              const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
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
				{           //bc
					scoring_per_pair.emplace_back(scoring_angle, current_pair);
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
            if(!b_robust_initialization_of_imu_)
            {
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
            
          }
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
      initial_score_log.open(stlplus::create_filespec(sOut_directory_,"initial_score.txt"));
      initial_score_log<<"view i and view j :score\n";
      for(const auto& pair: scoring_per_pair)
      {
          initial_score_log<<pair.second.first<<"-"<<pair.second.second<<":   "<<pair.first<<"\n";
      }
      initial_score_log.close();
  }
  
  /////////bc debug end///////
  if (!scoring_per_pair.empty())
  {
    initial_pair = scoring_per_pair.begin()->second;
    return true;
  }
  return false;
}

/// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
bool SequentialSfMReconstructionEngine_Robust_Initialization::RobustMakeInitialPair3D(const Pair & current_pair)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  std::cout<<"Robust Make InitialPair3D\n";
  const uint32_t
    I = std::min(current_pair.first, current_pair.second),
    J = std::max(current_pair.first, current_pair.second);

  if (sfm_data_.GetViews().count(I) == 0 ||
      sfm_data_.GetViews().count(J) == 0)
  {
    std::cout<<"missing view "<<I<<" or "<<J<<"\n";
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
    std::cout<<"missing intrinsics "<<I<<" or "<<J<<"\n";
    return false;
  }

  const auto
    * cam_I = iterIntrinsic_I->second.get(),
    * cam_J = iterIntrinsic_J->second.get();
  if (!cam_I || !cam_J)
  {
    std::cout<<"unvalid intrinsics "<<I<<" or "<<J<<"\n";
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
	  ransac_file << "Pair "<<I<<" "<<J<<"\n";
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
    std::ofstream ransac_bearinginput_file(stlplus::create_filespec(sOut_directory_, "ransac_bearinginput.txt"));
	  ransac_bearinginput_file << "Pair "<<I<<" "<<J<<"\n";
    ransac_bearinginput_file << "image size I(w,h):" << imageSize_I.first << " " << imageSize_I.second << "\n";
	  
	  ransac_bearinginput_file << "intrinsics I:\n";
	  
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  ransac_bearinginput_file << KI(i, j) << " ";
		  }
		  ransac_bearinginput_file << "\n";
	  }
	  ransac_bearinginput_file << "image size J(w,h):" << imageSize_J.first << " " << imageSize_J.second << "\n";
	  ransac_bearinginput_file << "intrinsics J:\n";
	  
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++)
		  {
			  ransac_bearinginput_file << KJ(i, j) << " ";
		  }
		  ransac_bearinginput_file << "\n";
	  }
    Mat3X bearing1 = (*cam_I)(xI);
    Mat3X bearing2 = (*cam_J)(xJ);
	  ransac_bearinginput_file << "bearing_xi bearing_yi bearing_xj bearing_yj:\n";
	  for (int i = 0; i < n; i++)
	  {
		  ransac_bearinginput_file << bearing1(0, i) << " " << bearing1(1, i) << " ";
		  ransac_bearinginput_file << bearing2(0, i) << " " << bearing2(1, i) << "\n";
	  }
	  ransac_bearinginput_file.close();
  }
  bool initialization_flag = false;
  if (!b_robust_initialization_of_imu_)
  {
  if (robustRelativePose_MultiModel(
    cam_I, cam_J, xI, xJ, relativePose_info, relativePose_essentialinfo,
    relativePose_homographyinfo,imageSize_I, imageSize_J, 4096))
      {
        initialization_flag = true;
      }
  }
  else
  {
      const Pose3 imu_pose_I = imu_data_.poses.at(I);
      const Pose3 imu_pose_J = imu_data_.poses.at(J);
      Pose3 imu_relative_pose = imu_pose_J * imu_pose_I.inverse();
      if(robustRelativePose_IMU(
        cam_I, cam_J,
        xI, xJ, relativePose_info,imu_relative_pose,
        imageSize_I, imageSize_J,
        4096))
      {
          initialization_flag = true;
      }
  }
  if (!initialization_flag)
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
  // Bound min precision at 1 pix.
  relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

  const bool bRefine_using_BA = true;
  if (bRefine_using_BA)
  {
    // Refine the defined scene
    SfM_Data tiny_scene;
    tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
    tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
    tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
    tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

    // Init poses
    const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

    // Init structure
    Landmarks & landmarks = tiny_scene.structure;

    for (const auto & track_iterator : map_tracksCommon)
    {
      // Get corresponding points
      auto iter = track_iterator.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second;

      const Vec2
        x1 = features_provider_->feats_per_view[I][i].coords().cast<double>(),
        x2 = features_provider_->feats_per_view[J][j].coords().cast<double>();

      Vec3 X;
      if (Triangulate2View(
            Pose_I.rotation(),
            Pose_I.translation(),
            (*cam_I)(cam_I->get_ud_pixel(x1)),
            Pose_J.rotation(),
            Pose_J.translation(),
            (*cam_J)(cam_J->get_ud_pixel(x2)),
            X,
            triangulation_method_))
      {
        Observations obs;
        obs[view_I->id_view] = Observation(x1, i);
        obs[view_J->id_view] = Observation(x2, j);
        landmarks[track_iterator.first].obs = std::move(obs);
        landmarks[track_iterator.first].X = X;
      }
    }
    Save(tiny_scene, stlplus::create_filespec(sOut_directory_, "initialPair.ply"), ESfM_Data(ALL));

    // - refine only Structure and Rotations & translations (keep intrinsic constant)
    Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.Adjust(tiny_scene,
        Optimize_Options
        (
          Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL) // Adjust structure
        )
      )
    {
      return false;
    }

    // Save computed data
    const Pose3 pose_I = sfm_data_.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
    const Pose3 pose_J = sfm_data_.poses[view_J->id_pose] = tiny_scene.poses[view_J->id_pose];
    map_ACThreshold_.insert({I, relativePose_info.found_residual_precision});
    map_ACThreshold_.insert({J, relativePose_info.found_residual_precision});
    set_remaining_view_id_.erase(view_I->id_view);
    set_remaining_view_id_.erase(view_J->id_view);

    // List inliers and save them
    std::ofstream initial_triangulation_log;////bc
    if(bdebug)
    {
		std::ofstream ba_initial_essential_log(stlplus::create_filespec(sOut_directory_, "ba_initial_essentials.txt"));
		ba_initial_essential_log << "ba relative pose:\n";
		Pose3 relative_poseji = pose_J * pose_I.inverse();
		Mat3 relative_rji = relative_poseji.rotation();
		Vec3 relative_tji = relative_poseji.translation();
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				ba_initial_essential_log << relative_rji(i, j) << " ";
			}
			ba_initial_essential_log << relative_tji(i) << "\n";
		}
		Mat3 tji_ni;
		tji_ni << 0,                -relative_tji(2), relative_tji(1),
			      relative_tji(2),  0,                -relative_tji(0),
			      -relative_tji(1), relative_tji(0),  0;
		Mat3 essential_matrix = tji_ni * relative_rji;
		ba_initial_essential_log << "ba essential matrix:\n";
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				ba_initial_essential_log << essential_matrix(i, j) << " ";
			}
			ba_initial_essential_log << "\n";
		}
		ba_initial_essential_log << "essential_error\n";
		std::vector<double> essential_error;
		for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
		{
			openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
			std::advance(iterT, inlier_idx);
			tracks::submapTrack::const_iterator iter = iterT->second.begin();
			assert(map_tracksCommon.at(inlier_idx).at(I) == iter->second);
			assert(map_tracksCommon.at(inlier_idx).at(J) == (++iter)->second);
			const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
			const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
			Vec2 norm_featI=cam_I->ima2cam(featI);
			Vec2 norm_featJ=cam_J->ima2cam(featJ);
			Vec3 norm_featI_3, norm_featJ_3;
			norm_featI_3 << norm_featI(0), norm_featI(1), 1.0;
			norm_featJ_3 << norm_featJ(0), norm_featJ(1), 1.0;
			Vec error = norm_featJ_3.transpose() * essential_matrix * norm_featI_3;
			essential_error.push_back(error(0));
			
		}
		std::sort(essential_error.begin(), essential_error.end());
		for (size_t i = 0; i < essential_error.size(); i++)
		{
			ba_initial_essential_log << essential_error[i] << "\n";
		}
		ba_initial_essential_log.close();
      initial_triangulation_log.open(stlplus::create_filespec(sOut_directory_,"initial_triangulation_status.txt"));
	  initial_triangulation_log << view_I->id_view << "-" << view_J->id_view << "\n";
	}
    for (const auto & landmark_entry : tiny_scene.GetLandmarks())
    {
      const IndexT trackId = landmark_entry.first;
      const Landmark & landmark = landmark_entry.second;
      const Observations & obs = landmark.obs;
      Observations::const_iterator
        iterObs_xI = obs.find(view_I->id_view),
        iterObs_xJ = obs.find(view_J->id_view);

      const Observation & ob_xI = iterObs_xI->second;
      const Observation & ob_xJ = iterObs_xJ->second;
      const Vec2
        ob_xI_ud = cam_I->get_ud_pixel(ob_xI.x),
        ob_xJ_ud = cam_J->get_ud_pixel(ob_xJ.x);

      const double angle = AngleBetweenRay(
        pose_I, cam_I, pose_J, cam_J, ob_xI_ud, ob_xJ_ud);
      const Vec2 residual_I = cam_I->residual(pose_I(landmark.X), ob_xI.x);
      const Vec2 residual_J = cam_J->residual(pose_J(landmark.X), ob_xJ.x);
      if (angle > 2.0 &&
          CheiralityTest((*cam_I)(ob_xI_ud), pose_I,
                         (*cam_J)(ob_xJ_ud), pose_J,
                         landmark.X) &&
          residual_I.norm() < relativePose_info.found_residual_precision &&
          residual_J.norm() < relativePose_info.found_residual_precision)
      {
        sfm_data_.structure[trackId] = landmarks[trackId];
		if (bdebug)
		{
			initial_triangulation_log << ob_xI.id_feat << " " << ob_xJ.id_feat 
				<<"("<< ob_xI.x(0)<<" "<< ob_xI.x(1)<<" "<< ob_xJ.x(0)<<" "<< ob_xJ.x(1)<<")"
				<< ": triangulated\n";  //bc
		}
      }
      else if(bdebug)  //bc
      {
        initial_triangulation_log << ob_xI.id_feat<<" "<<ob_xJ.id_feat<<": ";  //bc
        if(!(angle > 2.0))
        {
          initial_triangulation_log<<" angle <= 2.0 ";
        }
        else if(!CheiralityTest((*cam_I)(ob_xI_ud), pose_I,
                         (*cam_J)(ob_xJ_ud), pose_J,
                         landmark.X))
        {
          initial_triangulation_log<<" CheiralityTest failed ";
        }
        else if(!(residual_I.norm() < relativePose_info.found_residual_precision ))
        {
          initial_triangulation_log<<" Reprjection of "<<ob_xI.id_feat<<" is too large ";
        }
        else if(!(residual_J.norm() < relativePose_info.found_residual_precision))
        {
          initial_triangulation_log<<" Reprjection of "<<ob_xJ.id_feat<<" is too large ";
        }
        initial_triangulation_log<<"\n";
      }
      
    }
    initial_triangulation_log.close();
    // Save outlier residual information
    Histogram<double> histoResiduals;
    std::cout << "\n"
      << "=========================\n"
      << " MSE Residual InitialPair Inlier:\n";
    ComputeResidualsHistogram(&histoResiduals);
    std::cout << "=========================" << std::endl;

    if (!sLogging_file_.empty())
    {
      using namespace htmlDocument;
      html_doc_stream_->pushInfo(htmlMarkup("h1","Essential Matrix."));
      std::ostringstream os;
      os << std::endl
        << "-------------------------------" << "<br>"
        << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
        << view_I->s_Img_path << ","
        << view_J->s_Img_path << "<br>"
        << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
        << "-- Resection status: " << "OK" << "<br>"
        << "-- Nb points used for robust Essential matrix estimation: "
        << xI.cols() << "<br>"
        << "-- Nb points validated by robust estimation: "
        << sfm_data_.structure.size() << "<br>"
        << "-- % points validated: "
        << sfm_data_.structure.size()/static_cast<float>(xI.cols())
        << "<br>"
        << "-------------------------------" << "<br>";
      html_doc_stream_->pushInfo(os.str());

      html_doc_stream_->pushInfo(htmlMarkup("h2",
        "Residual of the robust estimation (Initial triangulation). Thresholded at: "
        + toString(relativePose_info.found_residual_precision)));

      html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of residuals"));

      const std::vector<double> xBin = histoResiduals.GetXbinsValue();
      const auto range = autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

      htmlDocument::JSXGraphWrapper jsxGraph;
      jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
      jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
      jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
        relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
      jsxGraph.UnsuspendUpdate();
      jsxGraph.setViewport(range);
      jsxGraph.close();
      html_doc_stream_->pushInfo(jsxGraph.toStr());

      html_doc_stream_->pushInfo("<hr>");

      std::ofstream htmlFileStream( std::string(stlplus::folder_append_separator(sOut_directory_) +
        "Reconstruction_Report.html").c_str());
      htmlFileStream << html_doc_stream_->getDoc();
    }
  }
  return !sfm_data_.structure.empty();
}

double SequentialSfMReconstructionEngine_Robust_Initialization::ComputeResidualsHistogram(Histogram<double> * histo)
{
  // Collect residuals for each observation
  std::vector<float> vec_residuals;
  vec_residuals.reserve(sfm_data_.structure.size());
  for (const auto & landmark_entry : sfm_data_.GetLandmarks())
  {
    const Observations & obs = landmark_entry.second.obs;
    for (const auto & observation : obs)
    {
      const View * view = sfm_data_.GetViews().find(observation.first)->second.get();
      const Pose3 pose = sfm_data_.GetPoseOrDie(view);
      const auto intrinsic = sfm_data_.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose(landmark_entry.second.X), observation.second.x);
      vec_residuals.emplace_back( std::abs(residual(0)) );
      vec_residuals.emplace_back( std::abs(residual(1)) );
    }
  }
  // Display statistics
  if (vec_residuals.size() > 1)
  {
    float dMin, dMax, dMean, dMedian;
    minMaxMeanMedian<float>(vec_residuals.cbegin(), vec_residuals.cend(),
                            dMin, dMax, dMean, dMedian);
    if (histo)  {
      *histo = Histogram<double>(dMin, dMax, 10);
      histo->Add(vec_residuals.cbegin(), vec_residuals.cend());
    }

    std::cout << std::endl << std::endl;
    std::cout << std::endl
      << "SequentialSfMReconstructionEngine_General::ComputeResidualsMSE." << "\n"
      << "\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size() << std::endl
      << "\t-- Residual min:\t" << dMin << std::endl
      << "\t-- Residual median:\t" << dMedian << std::endl
      << "\t-- Residual max:\t "  << dMax << std::endl
      << "\t-- Residual mean:\t " << dMean << std::endl;

    return dMean;
  }
  return -1.0;
}




} // namespace sfm
} // namespace openMVG
