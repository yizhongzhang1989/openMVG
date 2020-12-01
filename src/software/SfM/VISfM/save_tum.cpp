//
// Created by v-xinli1 on 10/26/2020.
//

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include <fstream>

using namespace openMVG;
using namespace openMVG::sfm;

void SaveTumPose(
        const std::string& save_path,
        const std::vector<double>& times,
        const std::vector<Eigen::Vector3d>& trans,
        const std::vector<Eigen::Quaterniond>& rotats)
{
    std::ofstream my_file;
    my_file.open( save_path );
    for( int index = 0;index < trans.size();++index  )
    {
        auto &tran = trans[index];
        auto &rotat = rotats[index];

        my_file.precision(20);
        my_file
        << times[index] << " ";
        my_file.precision(5);
        my_file
        << tran(0) << " " << tran(1) << " " << tran(2) << " "
        << rotat.x() << " " << rotat.y() << " " << rotat.z() <<" " <<  rotat.w() << std::endl;
    }
    my_file.close();
}

int main(int argc , char** argv)
{
    CmdLine cmd;

    std::string sSfM_Data_Filename_In, sInputTimes;
    std::string sOutputPose_Oout;
    std::string sSfM_IMU_FileType;
    bool simu_mode = false;
    bool save_visual = false;

    cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
    cmd.add(make_option('o', sOutputPose_Oout, "output_pose_file"));
    cmd.add(make_option('t', sInputTimes, "input_time_file"));
    cmd.add( make_option('I', sSfM_IMU_FileType,"imu_filetype"));

    cmd.add( make_switch('v', "visual_save"));
//    cmd.add(make_option('t', sOutputPoint_Oout, "output_point_file"));

    try
    {
        if (argc == 1)  throw std::string(" Please input sfm_data.bib dic ");
        cmd.process(argc, argv);
    }
    catch (const std::string& s)
    {
        std::cerr << "ERROR" << std::endl;
        return EXIT_FAILURE;
    }

    Mat3 Ric;
    Vec3 tic;
    if( sSfM_IMU_FileType == std::string("EuRoc") )
    {
        Ric << 0.0148655429818, -0.999880929698, 0.00414029679422,
                0.999557249008, 0.0149672133247, 0.025715529948,
                -0.0257744366974, 0.00375618835797, 0.999660727178;
        tic << -0.0216401454975, -0.064676986768, 0.00981073058949;
    }
    else if ( sSfM_IMU_FileType == std::string("Mate20Pro") )
    {
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
        Ric.setIdentity();
        tic <<.0, .0, .0;

        simu_mode = true;
    }
    else
    {
        std::cerr <<"sSfM_IMU_FileType error  " << sSfM_IMU_FileType << std::endl;
        return 1;
    }

    save_visual = cmd.used('v');

    if(save_visual)
    {
        Ric.setIdentity();
        tic.setZero();
    }

    SfM_Data sfm_data;
    if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
    {
        std::cerr << std::endl
                  << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<double> times_in;
    {

        std::ifstream fin(sInputTimes, std::ifstream::in);
        int i=1;
        while( !fin.eof() )
        {
            double time;
            fin >> time;

            times_in.emplace_back( time  ); // second * 1000
//            if( sSfM_IMU_FileType == std::string("Mate20Pro") )
//            {
//                times.emplace_back( static_cast<IndexT>(time * 1000) ); // second * 1000
//            }
//            else if( sSfM_IMU_FileType == std::string("EuRoc") )
//            {
//                time -= 1403715 * 1e12;
//                times.emplace_back( static_cast<IndexT>(time / 1e6) ); // second * 1000
//            }
        }
        fin.close();
    }
    std::map<IndexT, double> times_map;
    int index = 0;
    for( auto&Id_View:sfm_data.views )
    {
        auto Id = Id_View.second->id_pose;
        times_map.emplace( Id, times_in[index++] );
    }

    auto poses = sfm_data.GetPoses();
    std::vector<Eigen::Vector3d> trans;
    std::vector<Eigen::Quaterniond> rotats;
    std::vector<double> times;
    for( auto &pose:poses)
    {
        double time;
        if( !simu_mode )
            time = times_map.at(pose.first) / 1e9;
        else
            time = times_map.at(pose.first);

        Eigen::Matrix3d Rcw = pose.second.rotation();
        Eigen::Matrix3d Rwc = Rcw.transpose();

        Eigen::Vector3d tcw = pose.second.translation();
        Eigen::Vector3d twc = -Rwc*tcw;

        Eigen::Vector3d tci = -Ric.transpose() * tic;
        Eigen::Matrix3d Rwi = Rwc * Ric.transpose();
        Eigen::Vector3d twi = twc + Rwc * tci;

        Eigen::Quaterniond Qwi( Rwi );
        times.push_back(time);
        trans.push_back(twi);
        rotats.push_back(Qwi);
    }

    SaveTumPose( sOutputPose_Oout, times, trans, rotats );

    return 0;
}