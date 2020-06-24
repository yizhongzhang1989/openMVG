// This file is part of OpenMVG_IMU , a branch of OpenMVG
// Author: Bao Chong
// Date:2020/06


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
  const size_t initial_max_iteration_count,
  const std::string & sloggingFile)
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

  if (!InitLandmarkTracks())
    return false;
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  // Initial pair choice
  openMVG::system::Timer initial_timer;
  if (initial_pair_ == Pair(0,0))
  {
    if (!RobustAutomaticInitialPairChoice(initial_pair_))
    {
      bool flag = Process_KnownTracks();   //turn the level of initialization down into original initialization
      return flag;
    }
  }
  std::cout << std::endl << " Total Initialization took (s): " << initial_timer.elapsed() << std::endl;
  // Else a starting pair was already initialized before

  // Initial pair Essential Matrix and [R|t] estimation.
  if (!RobustMakeInitialPair3D(initial_pair_))
    return false;
  //END(Author: BC)===================================================
  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  size_t resectionGroupIndex = 0;
  std::vector<uint32_t> vec_possible_resection_indexes;
  
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (const auto & iter : vec_possible_resection_indexes)
    {
      bImageAdded |= Resection(iter);
      set_remaining_view_id_.erase(iter);
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
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  if (badTrackRejector(4.0, 0))
  {
    eraseUnstablePosesAndObservations(sfm_data_);
  }
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
  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  size_t resectionGroupIndex = 0;
  std::vector<uint32_t> vec_possible_resection_indexes;
  
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (const auto & iter : vec_possible_resection_indexes)
    {
      bImageAdded |= Resection(iter);
      set_remaining_view_id_.erase(iter);
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
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  if (badTrackRejector(4.0, 0))
  {
    eraseUnstablePosesAndObservations(sfm_data_);
  }
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
          ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
          // Robust estimation of the relative pose
          RelativePose_MultiInfo relativePose_info;
					
          bool initialization_flag = false;
          if (!b_robust_initialization_of_imu_)
          {
            //using multiple geometry model in initialization
          if (robustRelativePose_MultiModel(
                cam_I, cam_J,
                xI, xJ, relativePose_info,
                {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
                initial_max_iteration_count_)
              && relativePose_info.vec_inliers.size() > iMin_inliers_count)
              {
                initialization_flag = true;
              }
          }
          else
          {
              //using imu validation in initialization
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
                  initialization_flag = true;
              }
          }
          //END(Author: BC)===================================================
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
				{           
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
  ////START(Author: BC)++++++++++++++++++++++++++++++++++++++++++++++
  // c. Robust estimation of the relative pose
  RelativePose_MultiInfo relativePose_info;
  
  const std::pair<size_t, size_t>
    imageSize_I(cam_I->w(), cam_I->h()),
    imageSize_J(cam_J->w(), cam_J->h());

  bool initialization_flag = false;
  if (!b_robust_initialization_of_imu_)
  {
    //using multiple geometry model in initialization
  if (robustRelativePose_MultiModel(
    cam_I, cam_J, xI, xJ, relativePose_info,imageSize_I, imageSize_J, 4096))
      {
        initialization_flag = true;
      }
  }
  else
  {
      //using imu validation in initialization
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
  //END(Author: BC)===================================================
  if (!initialization_flag)
  {
    std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
      << std::endl;
    return false;
  }

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
      }
      
    }

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
