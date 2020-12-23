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

#include "VISfM_src/sequential_VISfM.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "../Window_SfM/xin_time.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;


void PrintExtric(const SfM_Data & sfm_data)
{
    std::cout << "sfm_data.IG_Ric = \n" << sfm_data.IG_Ric << std::endl;
    std::cout << "sfm_data.IG_tic = \n" << sfm_data.IG_tic << std::endl;
}

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
    std::string sSfM_IMU_FileType = "Mate20Pro";
    std::string sSfM_Stamp_Filename;
    bool b_simu_g = false;


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
    cmd.add( make_option('I', sSfM_IMU_FileType,"imu_filetype"));
    cmd.add( make_option('T',sSfM_Stamp_Filename,"timestamp_file"));

    cmd.add( make_switch('g', "simu_gravity") );

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

    // EuRoc
//        0.0148655429818, -0.999880929698, 0.00414029679422, -0.0216401454975,
//         0.999557249008, 0.0149672133247, 0.025715529948, -0.064676986768,
//        -0.0257744366974, 0.00375618835797, 0.999660727178, 0.00981073058949,
//         0.0, 0.0, 0.0, 1.0

    b_simu_g = cmd.used('g');

    Mat3 Ric;
    Vec3 tic;
    if( sSfM_IMU_FileType == std::string("EuRoc") )
    {
// EuRoc
        Ric << 0.0148655429818, -0.999880929698, 0.00414029679422,
            0.999557249008, 0.0149672133247, 0.025715529948,
            -0.0257744366974, 0.00375618835797, 0.999660727178;
        tic << -0.0216401454975, -0.064676986768, 0.00981073058949;
    }
    else if ( sSfM_IMU_FileType == std::string("Mate20Pro") )
    {
//        //Mate20Pro
//        Ric << -0.00268725, -0.99990988, -0.0131532,
//                -0.99995582, 0.00280539, -0.00897172,
//                0.00900781, 0.01312851, -0.99987324;
//        tic << 0.01903381, -0.02204486, 0.00402214;

//        tic <<  0.0291997, -0.0430028,   0.109315;
//        Ric <<
//                -0.00952889,   -0.999334,  -0.0352139,
//                -0.999732,  0.00877768,   0.0214261,
//                -0.0211028,   0.0354086,    -0.99915;

        //Mate20Pro
//        Mat3 Rci;
//        Vec3 tci;
//
//        Rci << -0.00080983, -0.99991118,  0.01330307,
//            -0.99981724,  0.0010637,   0.01908794,
//            -0.01910039, -0.01328518, -0.9997293;
//        tci << 0.02532891, 0.03078696, 0.080411;
//
//        Ric = Rci.transpose();
//        tic = -Ric * tci;
//
//
//        std::cout << "tic = " << tic.transpose() << std::endl;
//        std::cout << "Ric = \n" << Ric << std::endl;
//        tic << 0.00807,
//        -0.00301,
//        -0.00273;
//        tic <<   0.07, 0., 0.02;//0.0294058, -0.00618251,   0.0419311;
        tic << 0.0294058, -0.00618251,   0.0419311;
        Ric <<
                0.0031626,   -0.999771,  -0.0211788,
                -0.999949, -0.00336582,  0.00956673,
                -0.00963582,   0.0211475,    -0.99973;

    }
    else if( sSfM_IMU_FileType == std::string("Simu") )
    {
        Ric.setIdentity();
        tic << 0., 0., 0.;
//        tic.setZero();
    }
    else
    {
        std::cerr <<"sSfM_IMU_FileType error  " << sSfM_IMU_FileType << std::endl;
        return 1;
    }



    Vec3 G;


    if( sSfM_IMU_FileType == std::string("Simu") )
    {
        if(b_simu_g)
            G = Vec3(0.,0.,9.8);
        else
            G = Vec3(0., 0., 0.);
    }
    else G = Vec3(0.,0.,9.8107);
    VIstaticParm::G_ = G;
    sfMData.IG_Ric = Ric;
    sfMData.IG_tic = tic;

    if( sSfM_IMU_FileType == std::string("EuRoc")  )
    {
        VIstaticParm::acc_n = 0.08;
        VIstaticParm::acc_w = 0.00004;
        VIstaticParm::gyr_n = 0.004;
        VIstaticParm::gyr_w = 2.0e-6;
    }
    else if( sSfM_IMU_FileType == std::string("Mate20Pro") )
    {
        VIstaticParm::acc_n = 0.01;
        VIstaticParm::acc_w = 0.0002;
        VIstaticParm::gyr_n = 0.005;
        VIstaticParm::gyr_w = 4.0e-6;

//        VIstaticParm::acc_n = 0.08;
//        VIstaticParm::acc_w = 0.00004;
//        VIstaticParm::gyr_n = 0.004;
//        VIstaticParm::gyr_w = 2.0e-6;

//        VIstaticParm::acc_n = 1.3061437477214574e-02;
//        VIstaticParm::acc_w = 9.7230832140122278e-04;
//        VIstaticParm::gyr_n = 9.7933408260869451e-04;
//        VIstaticParm::gyr_w = 1.7270393511376410e-05;
    }
    else
    {
        VIstaticParm::acc_n = 0.08;
        VIstaticParm::acc_w = 0.00004;
        VIstaticParm::gyr_n = 0.004;
        VIstaticParm::gyr_w = 2.0e-6;

//        VIstaticParm::acc_n = 0;
//        VIstaticParm::acc_w = 0;
//        VIstaticParm::gyr_n = 0;
//        VIstaticParm::gyr_w = 0;

//        VIstaticParm::acc_n = 1.0e-2;
//        VIstaticParm::acc_w = 1.0e-2;
//        VIstaticParm::gyr_n = 1.0e-2;
//        VIstaticParm::gyr_w = 1.0e-2;
    }


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
    std::cout << "match load begin"<< std::endl;
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

    std::cout << "match load over" << std::endl;

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

    std::vector<double> times;
    {

        std::ifstream fin(sSfM_Stamp_Filename, std::ifstream::in);
        if( !fin.is_open() ) return -1;
        int i=1;
        while( !fin.eof() )
        {
            double time;
            fin >> time;
//            std::string line;
//            getline(fin, line);
//            if( line.empty() )
//            {
////                std::cerr << "WTF"<< line << "WTF" << std::endl;
//                break;
//            }
////            double time = std::stod(line);
            if( sSfM_IMU_FileType == std::string("Mate20Pro") || sSfM_IMU_FileType == std::string("Simu") )
            {
                times.emplace_back( static_cast<unsigned long long int>(time * 1000) ); // second * 1000
            }
            else if( sSfM_IMU_FileType == std::string("EuRoc") )
            {
                time -= 1403715 * 1e12;
                time = static_cast<uint64_t>(time / 1e6);
                times.emplace_back( time ); // second * 1000
            }
//            static_cast<IndexT>(time / 1e6)
//             double time = std::stod(line);

////            std::cout << i++ << " " ;
//            std::cout.precision(20);
//            std::cout << time << std::endl;
//            std::cout << times.back() << std::endl;
        }
        fin.close();
    }

    std::cout << "time file read over" << std::endl;

    if( sSfM_IMU_Filename.empty() )
    {
        std::cerr << "not input sSfM_IMU_Filename " << std::endl;
        return EXIT_FAILURE;
    }
    // TODO xinli check IMU image dt
    std::shared_ptr<IMU_Dataset> imu_dataset = std::make_shared<IMU_Dataset>(sSfM_IMU_Filename, sSfM_IMU_FileType);
    std::cout << "imu file read over" << std::endl;
    std::cout << "imu time size for second: " << imu_dataset->vec_times.size() / 200 << std::endl;




//    if( sSfM_IMU_FileType == std::string( "Mate20Pro" ) )
//        imu_dataset->corect_time( times.back() );

//    IndexT dt = ;
    if(sSfM_IMU_FileType == std::string( "Mate20Pro" ))
//        imu_dataset->corect_dt( -0.26 * 1000 );
        imu_dataset->corect_dt( -0.2 * 1000 );

//    imu_dataset->corect_dt( 0.05 * 1000 );

//    imu_dataset->corect_dt( 50);

    // IMU_Dataset::sum_st_ = -50.1366

    IMU_Dataset::sum_st_ = 0.;

//    timeshift cam0 to imu0: [s] (t_imu = t_cam + shift)n sfm_data_.views.size() =

//    imu_dataset->corect_dt( -2 );

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


    std::string output_result_file = stlplus::create_filespec(sOutDir, "td_ex", ".txt");
    visfmEngine.output_result_file_ = output_result_file;
    XinTime time;
    if(visfmEngine.VI_Init(  ))
    {

        Save(visfmEngine.Get_SfM_Data(),
             stlplus::create_filespec(sOutDir, "sfm_visual_init_data", ".bin"),
             ESfM_Data(ALL));

        Save(visfmEngine.Get_SfM_Data(),
             stlplus::create_filespec(sOutDir, "sfm_visual_init_cloud_and_poses", ".ply"),
             ESfM_Data(ALL));
//        Save(visfmEngine.Get_SfM_Data(),
//             "/home/xinli/work/data/VI_visual_init.bin",
//             ESfM_Data(ALL));
        if(!visfmEngine.VI_align())
        {
            std::cout << "VI sfm align fail" << std::endl;

            bool align_ok = false;
            for(int i=0;i<10;++i)
            {
                if(visfmEngine.VI_Init_again(  ))
                {
                    if(visfmEngine.VI_align())
                    {
                        align_ok = true;
                        break;
                    }
                }
            }
            if(!align_ok)
                return EXIT_FAILURE;
        }

        visfmEngine.output_log_file_ = stlplus::create_filespec(sOutDir, "log", ".txt");
        {
//            Save(visfmEngine.Get_SfM_Data(),
//                 "/home/xinli/work/data/VI_visualIMU_init.bin",
//                 ESfM_Data(ALL));
//            return EXIT_SUCCESS;
            if(visfmEngine.Process_window())
//            if(visfmEngine.Process())
            {

//                time.print_time();
//                Save(visfmEngine.Get_SfM_Data(),
//                     "/home/xinli/work/data/Allresutl.bin",
//                     ESfM_Data(ALL));

                std::cout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;

                std::cout << "init ex\n";
                PrintExtric(sfMData);
                std::cout << "after oti" << std::endl;
                PrintExtric(visfmEngine.Get_SfM_Data());

                std::cout << "IMU_Dataset::sum_st_ = " << IMU_Dataset::sum_st_ << std::endl;


                {
                    std::ofstream fout( output_result_file, std::ofstream::app );
                    fout.precision(3);
//                    fout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;
                    time.toc();
                    fout << "FULL RECONSTRUCTION " << time.time_str();
                    fout << "IMU_Dataset::sum_st_ = " << IMU_Dataset::sum_st_ << std::endl;

                    auto sfm_data_opti = visfmEngine.Get_SfM_Data();

                    fout << "sfm_data.IG_Ric = \n" << sfm_data_opti.IG_Ric << std::endl;
                    fout << "sfm_data.IG_tic = \n" << sfm_data_opti.IG_tic << std::endl;
                    fout.close();
                }

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
