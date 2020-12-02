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
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <iostream>
#include <iomanip>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

class AppendSequentialSfMReconstructionEngine : public SequentialSfMReconstructionEngine
{
public:
    AppendSequentialSfMReconstructionEngine(const SfM_Data & sfm_data, const std::string & soutDirectory, const std::string & sloggingFile)
    : SequentialSfMReconstructionEngine(sfm_data,soutDirectory,sloggingFile)
    {
        ;
    }
    bool Process() override
    {
        using namespace openMVG::tracks;
        using namespace openMVG::features;
        //-------------------
        //-- Incremental reconstruction
        //-------------------

        if (!InitLandmarkTracks())
            return false;

        // reset landmarks
        Landmarks & structure = sfm_data_.structure;
        structure.clear();

        // Fill sfm_data with the computed tracks (no 3D yet)
        STLMAPTracks & map_selectedTracks = map_tracks_;
        IndexT idx(0);
        for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
             itTracks != map_selectedTracks.end();
             ++itTracks, ++idx)
        {
            const submapTrack & track = itTracks->second;
            structure[idx] = Landmark();
            Observations & obs = structure.at(idx).obs;
            for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
            {
                const size_t imaIndex = it->first;
                const size_t featIndex = it->second;
                const PointFeature & pt = features_provider_->feats_per_view.at(imaIndex)[featIndex];
                obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
            }
        }

        // Compute 3D position of the landmark of the structure by triangulation of the observations
        {
            openMVG::system::Timer timer;

            const IndexT trackCountBefore = sfm_data_.GetLandmarks().size();
            SfM_Data_Structure_Computation_Blind structure_estimator(true);
            structure_estimator.triangulate(sfm_data_);

            std::cout << "\n#removed tracks (invalid triangulation): " <<
                      trackCountBefore - IndexT(sfm_data_.GetLandmarks().size()) << std::endl;
            std::cout << std::endl << "  Triangulation took (s): " << timer.elapsed() << std::endl;

            // Export initial structure
            if (!sLogging_file_.empty())
            {
                Save(sfm_data_,
                     stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "initial_structure", "ply"),
                     ESfM_Data(EXTRINSICS | STRUCTURE));
            }
        }

        // bundle adjustment
        do
        {
            BundleAdjustment();
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
//        std::cout << "\nHistogram of residuals:\n" << h.ToString() << std::endl;

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

    double ComputeResidualsHistogram(Histogram<double> * histo)
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
};

int main(int argc, char **argv)
{
    using namespace std;
    std::cout << "Sequential/Incremental reconstruction" << std::endl
              << " Perform incremental SfM (Initial Pair Essential + Resection)." << std::endl
              << std::endl;

    CmdLine cmd;

    std::string sSfM_Data_Filename;
    std::string sMatchesDir, sMatchFilename;
    std::string sOutDir = "";

    std::string sIntrinsic_refinement_options = "ADJUST_ALL";
    int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
    bool b_use_motion_priors = false;
    int triangulation_method = static_cast<int>(ETriangulationMethod::DEFAULT);

    cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
    cmd.add( make_option('m', sMatchesDir, "matchdir") );
    cmd.add( make_option('M', sMatchFilename, "match_file") );
    cmd.add( make_option('o', sOutDir, "outdir") );
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

    if ( !isValid(static_cast<ETriangulationMethod>(triangulation_method)))
    {
        std::cerr << "\n Invalid triangulation method" << std::endl;
        return EXIT_FAILURE;
    }

    if ( !isValid(openMVG::cameras::EINTRINSIC(i_User_camera_model)) )
    {
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
    SfM_Data sfm_data;
    if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL)))
    {
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
    if (!feats_provider->load(sfm_data, sMatchesDir, regions_type))
    {
        std::cerr << std::endl
                  << "Invalid features." << std::endl;
        return EXIT_FAILURE;
    }
    // Matches reading
    std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
    if // Try to read the provided match filename or the default one (matches.f.txt/bin)
    (
        !(matches_provider->load(sfm_data, sMatchFilename) ||
          matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.f.txt")) ||
          matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.f.bin")))
    )
    {
        std::cerr << std::endl
                  << "Invalid matches file." << std::endl;
        return EXIT_FAILURE;
    }

    if (sOutDir.empty())
    {
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
    // Sequential reconstruction process
    //---------------------------------------

    openMVG::system::Timer timer;
    AppendSequentialSfMReconstructionEngine sfmEngine(
            sfm_data,
            sOutDir,
            stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

    // Configure the features_provider & the matches_provider
    sfmEngine.SetFeaturesProvider(feats_provider.get());
    sfmEngine.SetMatchesProvider(matches_provider.get());

    // Configure reconstruction parameters
    sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
    sfmEngine.SetUnknownCameraType(EINTRINSIC(i_User_camera_model));
    b_use_motion_priors = cmd.used('P');
    sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);
    sfmEngine.SetTriangulationMethod(static_cast<ETriangulationMethod>(triangulation_method));

    if (sfmEngine.Process())
    {

        std::cout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;

        std::cout << "...Generating SfM_Report.html" << std::endl;
        Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
                            stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

        //-- Export to disk computed scene (data & visualizable results)
        std::cout << "...Export SfM_Data to disk." << std::endl;
        Save(sfmEngine.Get_SfM_Data(),
             stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
             ESfM_Data(ALL));

        Save(sfmEngine.Get_SfM_Data(),
             stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
             ESfM_Data(ALL));

        return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}
