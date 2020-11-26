//
// Created by v-xinli1 on 11/23/2020.
//

#include "window_sequential_SfM.hpp"

#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
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

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#include <ceres/types.h>
#include <functional>
#include <iostream>
#include <utility>


namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;

    WindowSequentialSfMReconstructionEngine::WindowSequentialSfMReconstructionEngine(
            const SfM_Data & sfm_data,
            const std::string & soutDirectory,
            const std::string & sloggingFile)
            :
            SequentialSfMReconstructionEngine(sfm_data, soutDirectory, sloggingFile)
    {

    }

    WindowSequentialSfMReconstructionEngine::~WindowSequentialSfMReconstructionEngine()
    {
        if (!sLogging_file_.empty())
        {
            // Save the reconstruction Log
            std::ofstream htmlFileStream(sLogging_file_.c_str());
            htmlFileStream << html_doc_stream_->getDoc();
        }
    }

    bool WindowSequentialSfMReconstructionEngine::Process()
    {
        //-------------------
        //-- Incremental reconstruction
        //-------------------

        if (!InitLandmarkTracks())
            return false;

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

        std::cout << "---------------------------------------\n"
                  << "initial_pair_.first = " << initial_pair_.first << "\n"
                  <<  "initial_pair_.second = " << initial_pair_.second << "\n"
                  << "---------------------------------------" << std::endl;

        // Initial pair Essential Matrix and [R|t] estimation.
        if (!MakeInitialPair3D(initial_pair_))
            return false;

        //Save(sfm_data_, stlplus::create_filespec(sOut_directory_, "initialPair", ".bin"), ESfM_Data(ALL));

        // Compute robust Resection of remaining images
        // - group of images will be selected and resection + scene completion will be tried
        size_t resectionGroupIndex = 0;
        std::vector<uint32_t> vec_possible_resection_indexes;
        while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
        {
            bool bImageAdded = false;
            // Add images to the 3D reconstruction
            sfm_data_.add_viewId_cur.clear();
            for (const auto & iter : vec_possible_resection_indexes)
            {
                bImageAdded |= Resection(iter);
                set_remaining_view_id_.erase(iter);
                sfm_data_.add_viewId_cur.insert(iter);
            }

            if (bImageAdded)
            {
                // Scene logging as ply for visual debug
                std::ostringstream os;
                os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
                Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

                //Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".bin"), ESfM_Data(ALL));


                if( resectionGroupIndex % 10 == 0 )
                {
                    std::cout << "start full BA " << std::endl;
                    do
                    {
                        BundleAdjustment( );
                    }
                    while (badTrackRejector(4.0, 50));
                }
                else
                {
                    std::cout << "start local BA " << std::endl;

                    SfM_Data local_scene;
                    std::set< IndexT > local_scene_viewId;
                    CreateLocalScene(local_scene, local_scene_viewId, sfm_data_.add_viewId_cur);

                    // Perform BA until all point are under the given precision
                    do
                    {
                        BundleAdjustmentWindows(local_scene);
                    }
                    while (badTrackRejectorWindow(local_scene,4.0, 50));

					for (auto& view_index : local_scene_viewId)
					{
						const View* view_I = sfm_data_.GetViews().at(view_index).get();
						sfm_data_.intrinsics[view_I->id_intrinsic] = local_scene.intrinsics[view_I->id_intrinsic];
						sfm_data_.poses[view_I->id_pose] = local_scene.poses[view_I->id_pose];
					}
					// update camera 3d points and 2d points
					for (auto& structure_landmark_it : local_scene.structure) {
						sfm_data_.structure[structure_landmark_it.first].X =
							local_scene.structure[structure_landmark_it.first].X;
					}
                }
                eraseUnstablePosesAndObservations(sfm_data_);

            }
            ++resectionGroupIndex;
        }

        do
        {
            BundleAdjustment( );
        }
        while (badTrackRejector(4.0, 50));
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

	bool WindowSequentialSfMReconstructionEngine::CreateLocalScene(
        SfM_Data&               local_scene, 
        std::set<IndexT>&       local_scene_viewId, 
        const std::set<IndexT>& add_viewId_cur
    ) {
        //  collect all views that is close to current views		
		for (auto& viewId : add_viewId_cur){
			{
				auto pose_it = sfm_data_.poses.find(viewId);
				int i = 0;
				while (pose_it != sfm_data_.poses.begin())
				{
					local_scene_viewId.insert(pose_it->first);
					pose_it--;
					if (i++ == 15) break;
				}
				if (i == 16) local_scene_viewId.insert(pose_it->first);
			}
			{
				auto pose_it = sfm_data_.poses.find(viewId);
				int i = 0;
				while (pose_it != sfm_data_.poses.end())
				{
					local_scene_viewId.insert(pose_it->first);
					pose_it++;
					if (i++ == 15) break;
				}
			}
		}

        //  add views to local_scene
		for (auto viewId : local_scene_viewId) {
			const View* view_I = sfm_data_.GetViews().at(viewId).get();
			local_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
			local_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
			local_scene.poses.insert(*sfm_data_.GetPoses().find(view_I->id_pose));
			assert(view_I->id_view == view_I->id_pose);
		}

		//std::cout << "local_scene_viewId.size() = " << local_scene_viewId.size() << std::endl;

		for (auto& structure_landmark_it : sfm_data_.structure) {
			const Observations& obs = structure_landmark_it.second.obs;
			Observations local_obs;
			bool f_use_landmark = false;
			for (const auto& obs_it : obs) {
				const View* view = sfm_data_.views.at(obs_it.first).get();
				if (local_scene_viewId.count(view->id_view) == 1)
				{
					local_obs[obs_it.first] = obs_it.second;
					f_use_landmark = true;
				}
			}
			if (f_use_landmark) {
				local_scene.structure[structure_landmark_it.first].obs = std::move(local_obs);
				local_scene.structure[structure_landmark_it.first].X = sfm_data_.structure[structure_landmark_it.first].X;
			}
		}

        return true;
	}

    bool WindowSequentialSfMReconstructionEngine::BundleAdjustmentWindows(SfM_Data &sfm_data)
    {
        Bundle_Adjustment_Ceres::BA_Ceres_options options;
        if ( sfm_data_.GetPoses().size() > 100 &&
             (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
                )
            // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
        {
            options.preconditioner_type_ = ceres::JACOBI;
            options.linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
        else
        {
            options.linear_solver_type_ = ceres::DENSE_SCHUR;
        }
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
                ( ReconstructionEngine::intrinsic_refinement_options_,
                  Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
                  Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
                  Control_Point_Parameter(),
                  this->b_use_motion_prior_
                );
        return bundle_adjustment_obj.Adjust(sfm_data, ba_refine_options);
    }

    bool WindowSequentialSfMReconstructionEngine::badTrackRejectorWindow(SfM_Data &sfm_data, double dPrecision, size_t count)
    {
        const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(sfm_data, dPrecision, 2);
        const size_t nbOutliers_angleErr = 0;//RemoveOutliers_AngleError(sfm_data, 2.0);


        std::cout << "nbOutliers_residualErr = " << nbOutliers_residualErr << "\n"
                  << "dPrecision = " << dPrecision << "\n"
                  << "nbOutliers_angleErr = " << nbOutliers_angleErr << "\n"
                  << "count = " << count << std::endl;
        return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
    }

    double WindowSequentialSfMReconstructionEngine::ComputeResidualsHistogram(Histogram<double> *histo)
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
                      << "SequentialSfMReconstructionEngine::ComputeResidualsMSE." << "\n"
                      << "\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size() << std::endl
                      << "\t-- Residual min:\t" << dMin << std::endl
                      << "\t-- Residual median:\t" << dMedian << std::endl
                      << "\t-- Residual max:\t "  << dMax << std::endl
                      << "\t-- Residual mean:\t " << dMean << std::endl;

            return dMean;
        }
        return -1.0;
    }

}
}