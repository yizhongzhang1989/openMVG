//
// Created by root on 10/12/20.
//

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/multiview/solver_resection.hpp"
#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"

#include "software/SfM/addIMU/imu_integrator/sfm_imu.hpp"
#include "VISfM_src/sequential_VISfM.hpp"

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
    std::pair<std::string,std::string> initialPairString("","");
    std::string sIntrinsic_refinement_options = "ADJUST_ALL";
    int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
    bool b_use_motion_priors = false;
    int triangulation_method = static_cast<int>(ETriangulationMethod::DEFAULT);

    int selection_method = static_cast<int>(resection::SelectionMethod::DEFAULT);
    std::string sSfM_IMU_Filename;
    std::string sSfM_Stamp_Filename;


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

    cmd.add(make_option('s', selection_method, "selection_method"));
    cmd.add( make_option('u', sSfM_IMU_Filename,"imu_file"));
    cmd.add( make_option('T',sSfM_Stamp_Filename,"timestamp_file"));

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
    SfM_Data sfMData;
    if (!Load(sfMData, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
        std::cerr << std::endl
                  << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
    }

    Mat3 Ric;
    Ric << -0.00268725, -0.99990988, -0.0131532,
            -0.99995582, 0.00280539, -0.00897172,
            0.00900781, 0.01312851, -0.99987324;
    Vec3 tic( 0.01903381, -0.02204486, 0.00402214 );
    Vec3 G(0.,0.,9.8107);
    VIstaticParm::G_ = G;
    sfMData.IG_Ric = Ric;
    sfMData.IG_tic = tic;

//    // Load input imu_Data
//    SfM_IMU imu_data;
//    if(!sSfM_IMU_Filename.empty() && stlplus::file_exists(sSfM_IMU_Filename))
//    {
//        imu_data.ReadFromCSV(sSfM_IMU_Filename);
//    }

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
    if (!feats_provider->load(sfMData, sMatchesDir, regions_type)) {
        std::cerr << std::endl
                  << "Invalid features." << std::endl;
        return EXIT_FAILURE;
    }
    // Matches reading
    std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
    if // Try to read the provided match filename or the default one (matches.f.txt/bin)
            (
            !(matches_provider->load(sfMData, sMatchFilename) ||
              matches_provider->load(sfMData, stlplus::create_filespec(sMatchesDir, "matches.f.txt")) ||
              matches_provider->load(sfMData, stlplus::create_filespec(sMatchesDir, "matches.f.bin")))
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

    std::vector<IndexT> times;
    {

        std::ifstream fin(sSfM_Stamp_Filename, std::ifstream::in);
        int i=1;
        while( !fin.eof() )
        {
            std::string line;
            getline(fin, line);
            if( line.empty() )
            {
//                std::cerr << "WTF"<< line << "WTF" << std::endl;
                break;
            }
            double time = std::stod(line);
            times.emplace_back( static_cast<IndexT>(time * 1000) ); // second * 1000
//            std::cout << i++ << " " ;
//            std::cout << time << std::endl;
        }
        fin.close();
    }

    if( sSfM_IMU_Filename.empty() )
    {
        std::cerr << "not input sSfM_IMU_Filename " << std::endl;
        return EXIT_FAILURE;
    }
    // TODO xinli check IMU image dt
    std::shared_ptr<IMU_Dataset> imu_dataset = std::make_shared<IMU_Dataset>(sSfM_IMU_Filename);
//    imu_dataset->corect_time( times.back() );
//    imu_dataset->corect_dt( 2 );

    //---------------------------------------
    // Sequential reconstruction process
    //---------------------------------------

    openMVG::system::Timer timer;
    SequentialVISfMReconstructionEngine visfmEngine(
            sfMData,
            sOutDir,
            stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

    // Configure the features_provider & the matches_provider
    visfmEngine.SetFeaturesProvider(feats_provider.get());
    visfmEngine.SetMatchesProvider(matches_provider.get());
    visfmEngine.SetTimeStamp(times);
    visfmEngine.SetIMUDataset(imu_dataset);

    // Configure reconstruction parameters
    visfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
    visfmEngine.SetUnknownCameraType(EINTRINSIC(i_User_camera_model));
    b_use_motion_priors = cmd.used('P');
    visfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);
    visfmEngine.SetTriangulationMethod(static_cast<ETriangulationMethod>(triangulation_method));


    // Handle Initial pair parameter
    if (!initialPairString.first.empty() && !initialPairString.second.empty())
    {
        Pair initialPairIndex;
        if (!computeIndexFromImageNames(sfMData, initialPairString, initialPairIndex))
        {
            std::cerr << "Could not find the initial pairs <" << initialPairString.first
                      <<  ", " << initialPairString.second << ">!\n";
            return EXIT_FAILURE;
        }
        visfmEngine.setInitialPair(initialPairIndex);
    }


    if(visfmEngine.VI_Init(  ))
    {

        if(!visfmEngine.VI_align())
        {
            std::cerr << "VI sfm align fail" << std::endl;
            return EXIT_FAILURE;
        }
        else{
            Save(visfmEngine.Get_SfM_Data(),
                 "/home/xinli/work/data/VI_visual_init.bin",
                 ESfM_Data(ALL));
            return EXIT_SUCCESS;
            if(visfmEngine.Process())
            {
                std::cout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;

                std::cout << "...Generating SfM_Report.html" << std::endl;
                Generate_SfM_Report(visfmEngine.Get_SfM_Data(),
                                    stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

                //-- Export to disk computed scene (data & visualizable results)
                std::cout << "...Export SfM_Data to disk." << std::endl;
                Save(visfmEngine.Get_SfM_Data(),
                     stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
                     ESfM_Data(ALL));

                Save(visfmEngine.Get_SfM_Data(),
                     stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
                     ESfM_Data(ALL));

                return EXIT_SUCCESS;
            }
//            return EXIT_SUCCESS;
        }
    }
    else
    {
        std::cerr << "VI sfm init fail" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_FAILURE;
}
