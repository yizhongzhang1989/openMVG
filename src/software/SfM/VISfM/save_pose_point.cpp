//
// Created by root on 10/21/20.
//

//
// Created by root on 9/28/20.
//

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include <fstream>

using namespace openMVG;
using namespace openMVG::sfm;

void SavePose( const std::string& save_path,
               const std::vector<Eigen::Vector3d>& trans,
               const std::vector<Eigen::Quaterniond>& rotats)
{
    std::ofstream my_file;
    my_file.open( save_path );
    for( int index = 0;index < trans.size();++index  )
    {
        auto &tran = trans[index];
        auto &rotat = rotats[index];

        my_file << tran(0) << " " << tran(1) << " " << tran(2) << " "
                << rotat.w() << " " << rotat.x() << " " << rotat.y() << " " << rotat.z() << std::endl;
    }
    my_file.close();
}

void SavePoints( const std::string& save_path,
                 const std::vector<Eigen::Vector3d>& Points)
{
    std::ofstream my_file;
    my_file.open( save_path );
    for( auto& Point:Points  )
    {
        my_file << Point(0) << " " << Point(1) << " " << Point(2) << std::endl;
    }
    my_file.close();
}

int main(int argc , char** argv)
{
    CmdLine cmd;

    std::string sSfM_Data_Filename_In;
    std::string sOutputPose_Oout, sOutputPoint_Oout;

    cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
    cmd.add(make_option('c', sOutputPose_Oout, "output_pose_file"));
    cmd.add(make_option('p', sOutputPoint_Oout, "output_point_file"));

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

    SfM_Data sfm_data;
    if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
    {
        std::cerr << std::endl
                  << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
    }

    auto poses = sfm_data.GetPoses();
    std::vector<Eigen::Vector3d> trans;
    std::vector<Eigen::Quaterniond> rotats;
    for( auto &pose:poses)
    {
        pose.first;
        Eigen::Quaterniond Qwc( pose.second.rotation() );
        Eigen::Vector3d twc( pose.second.translation() );
        trans.push_back(twc);
        rotats.push_back(Qwc);
    }

    SavePose( sOutputPose_Oout, trans, rotats );

    std::vector<Eigen::Vector3d> Points;
    auto Landmarks = sfm_data.GetLandmarks();
    for(auto &landmark:Landmarks)
    {
        Points.push_back( landmark.second.X );
    }
    SavePoints( sOutputPoint_Oout, Points );

    return 0;
}