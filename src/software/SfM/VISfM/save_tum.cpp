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
    bool simu_mode = false;

    cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
    cmd.add(make_option('o', sOutputPose_Oout, "output_pose_file"));
    cmd.add(make_option('t', sInputTimes, "input_time_file"));
    cmd.add( make_switch('s', "simu_mode"));
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


    simu_mode = cmd.used('s');

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

        Eigen::Quaterniond Qwc( Rwc );
        times.push_back(time);
        trans.push_back(twc);
        rotats.push_back(Qwc);
    }

    SaveTumPose( sOutputPose_Oout, times, trans, rotats );

    return 0;
}