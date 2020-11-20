//
// Created by v-xinli1 on 11/19/2020.
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
#include <random>


using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;


bool ReadPoses( const std::string& sSfM_Pose_GT, Poses& poses )
{
    std::ifstream fin(sSfM_Pose_GT, std::ifstream::in);
    if( !fin.is_open() ) return false;
    int index = 0;
    while( !fin.eof() )
    {
        double data[8];
        for(double & i : data) fin >> i;
        Eigen::Vector3d twc( data[1], data[2], data[3] );
        Eigen::Quaterniond Qwc( data[7], data[4], data[5], data[6] );
        Mat3 rotation = Qwc.toRotationMatrix().transpose();
        Pose3 pose( rotation, twc );
        poses.insert( std::make_pair( index, pose ) );
//        std::cout << index << " " <<  twc.transpose() << std::endl;
        index++;
        if( index == 200 ) break;
    }
    fin.close();
}

int main(int argc, char **argv)
{
    CmdLine cmd;

    std::string sSfM_Data_Filename;
    std::string sOutDir = "";

    std::string sSfM_IMU_Filename;
    std::string sSfM_IMU_FileType = "Mate20Pro";
    std::string sSfM_Stamp_Filename;
    std::string sSfM_Pose_GT;


    cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
    cmd.add( make_option('o', sOutDir, "outdir") );
    cmd.add( make_option('p', sSfM_Pose_GT, "poseGT_file") );

    cmd.add( make_option('u', sSfM_IMU_Filename,"imu_file"));
    cmd.add( make_option('I', sSfM_IMU_FileType,"imu_filetype"));
    cmd.add( make_option('T',sSfM_Stamp_Filename,"timestamp_file"));


    try {
        if (argc == 1) throw std::string("Invalid parameter.");
        cmd.process(argc, argv);
    } catch (const std::string& s) {
        std::cerr << s << std::endl;
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

        //Mate20Pro
        Mat3 Rci;
        Vec3 tci;

        Rci << -0.00080983, -0.99991118,  0.01330307,
                -0.99981724,  0.0010637,   0.01908794,
                -0.01910039, -0.01328518, -0.9997293;
        tci << 0.02532891, 0.03078696, 0.080411;

        Ric = Rci.transpose();
        tic = -Ric * tci;

    }
    else if( sSfM_IMU_FileType == std::string("Simu") )
    {
//        Mat3 Rci;
//        Vec3 tci;
//        Rci << -0.00080983, -0.99991118,  0.01330307,
//                -0.99981724,  0.0010637,   0.01908794,
//                -0.01910039, -0.01328518, -0.9997293;
//        tci << 0.02532891, 0.03078696, 0.080411;
//        Ric = Rci.transpose();
//        tic = -Ric * tci;

//        Eigen::Matrix4d Tcc = Eigen::Matrix4d::Identity();
//        Eigen::Vector3d tcc = Eigen::Vector3d(nt_(generator_), nt_(generator_), nt_(generator_));



//        Eigen::Matrix3d Rcc = Eigen::Matrix3d::Identity();
//        Rcc = Eigen::AngleAxisd( ( 3 )/180.*M_PI , Eigen::Vector3d::UnitZ())
//              * Eigen::AngleAxisd( ( 3 )/180.*M_PI , Eigen::Vector3d::UnitY())
//              * Eigen::AngleAxisd( ( 3 )/180.*M_PI , Eigen::Vector3d::UnitX());

        Ric.setIdentity();

//        Ric = Ric * Rcc;
        tic <<.0, .0, .0;
//        tic <<.01, .01, .01;
//        tic <<.01, .01, .01;
//        tic <<.05,.05,.05;
    }
    else
    {
        std::cerr <<"sSfM_IMU_FileType error  " << sSfM_IMU_FileType << std::endl;
        return 1;
    }

    Vec3 G;


    if( sSfM_IMU_FileType == std::string("Simu") ) G = Vec3(0.,0.,0);
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
        VIstaticParm::acc_n = 1.3061437477214574e-02;
        VIstaticParm::acc_w = 9.7230832140122278e-04;
        VIstaticParm::gyr_n = 9.7933408260869451e-04;
        VIstaticParm::gyr_w = 1.7270393511376410e-05;
    }
    else
    {
        VIstaticParm::acc_n = 0.0;
        VIstaticParm::acc_w = 0.0;
        VIstaticParm::gyr_n = 0.0;
        VIstaticParm::gyr_w = 0.0;
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
        }
        fin.close();
    }
    std::cout << "times load end"<< std::endl;

    Poses poses;
    if( ReadPoses( sSfM_Pose_GT, poses ) )
    {
        std::cout << "poses.size() = " << poses.size() << std::endl;
    }
    else
    {
        std::cerr << "read pose GT error " << std::endl;
        return EXIT_FAILURE;
    }


    if( sSfM_IMU_Filename.empty() )
    {
        std::cerr << "not input sSfM_IMU_Filename " << std::endl;
        return EXIT_FAILURE;
    }
    std::shared_ptr<IMU_Dataset> imu_dataset = std::make_shared<IMU_Dataset>(sSfM_IMU_Filename, sSfM_IMU_FileType);
    std::cout << "imu file read over" << std::endl;
    std::cout << "imu time size for second: " << imu_dataset->vec_times.size() / 200 << std::endl;

    SequentialVISfMReconstructionEngine visfmEngine(
            sfMData,
            sOutDir,
            stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

    // Configure the features_provider & the matches_provider
    visfmEngine.SetTimeStamp(times);
    visfmEngine.SetIMUDataset(imu_dataset);
    visfmEngine.SetPoseGT(poses);


    if(visfmEngine.TestImuFactor())
    {
        std::cout << "...Export SfM_Data to disk." << std::endl;
        Save(visfmEngine.Get_SfM_Data(),
             stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
             ESfM_Data(ALL));

        return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;

}