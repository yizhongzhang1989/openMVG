#include "SimulationGenerator.h"
#include <unistd.h>
//#include <sys/types.h>
//#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>

// see https://www.cnblogs.com/renyu310/p/6485066.html for instructions of directory operation

namespace generator
{

bool SimulationGenerator::Generate(Simulation_Data& sfm_data, const SimulationConfig& cfg)
{
    // generate points
    MapPoints& map_points = sfm_data.map_points;
    unsigned int pt_id = 0;
    for(int i=0;i<cfg.n_points;i++)
    {
        map_points[pt_id] = MapPoint(pt_id,mpPointGenerator->Generate());
        pt_id++;
    }

    // generate poses
    KeyFrames& key_frames = sfm_data.key_frames;
    unsigned int kf_id = 0;
    for(int i=0;i<cfg.n_poses;i++)
    {
        key_frames[kf_id] = KeyFrame(kf_id,mpPoseGenerator->Generate());
        kf_id++;
    }
    std::cout<<"generation finished."<<std::endl;

    // construct relationships
    for(auto & key_frame : key_frames)
    {
        Pose& p = key_frame.second.pose;
//        Eigen::Matrix3d R = p.q.toRotationMatrix().transpose();
//        Eigen::Vector3d t = -R * p.t;
        for(auto & map_point : map_points)
        {
            Eigen::Vector3d& X = map_point.second.X;
//            Eigen::Vector3d Xc = R * X + t;
            Eigen::Vector3d Xc = p.q * X + p.t;
            Eigen::Vector2d xp = mpCamera->Project(Xc);
            if(mpCamera->isInView(xp))
            {
                key_frame.second.addObservation(map_point.second.Id);
                map_point.second.addObservation(key_frame.second.Id,xp);
            }
        }
    }
    std::cout<<"association finished."<<std::endl;
}

void SimulationGenerator::AddNoise(const Simulation_Data& sfm_data, Simulation_Data& sfm_data_noisy, NoiseConfig& cfg_noise)
{
    const MapPoints& map_points = sfm_data.map_points;
    const KeyFrames& key_frames = sfm_data.key_frames;
    MapPoints& map_points_noisy = sfm_data_noisy.map_points;
    KeyFrames& key_frames_noisy = sfm_data_noisy.key_frames;
    map_points_noisy.clear();
    key_frames_noisy.clear();

    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::normal_distribution<double> d_trans(0.0,cfg_noise.stddev_translation);
    std::normal_distribution<double> d_rot(0.0,cfg_noise.stddev_rotation);
    std::normal_distribution<double> d_point2d(0.0,cfg_noise.stddev_point2d);
    std::normal_distribution<double> d_point3d(0.0,cfg_noise.stddev_point3d);

    for(const auto & key_frame : key_frames)
    {
        Eigen::Quaterniond dq;
        if(cfg_noise.perturb_rotation)
        {
            dq = Eigen::AngleAxisd(d_rot(gen),Eigen::Vector3d::UnitZ())
                 *Eigen::AngleAxisd(d_rot(gen),Eigen::Vector3d::UnitY())
                 *Eigen::AngleAxisd(d_rot(gen),Eigen::Vector3d::UnitX());
        }
        else
        {
            dq.setIdentity();
        }

        Eigen::Vector3d dt;
        if(cfg_noise.perturb_translation)
        {
            dt = Eigen::Vector3d(d_trans(gen),d_trans(gen),d_trans(gen));
        }
        else
        {
            dt.setZero();
        }

        const Pose& pose = key_frame.second.pose;
        unsigned int id = key_frame.second.Id;
        Pose pose_noisy;
        pose_noisy.q = dq * pose.q;
        pose_noisy.t = dt + pose.t;
        key_frames_noisy[id] = KeyFrame(id,pose_noisy);
    }
    std::cout<<"add noise for pose finished."<<std::endl;

    for(const auto & map_point : map_points)
    {
        unsigned int id = map_point.second.Id;
//        const std::map<unsigned int, Eigen::Vector2d>& obs = map_point.second.obs;
        const STLMap<unsigned int, Eigen::Vector2d>& obs = map_point.second.obs;

        Eigen::Vector3d dX;
        if(cfg_noise.perturb_point3d)
        {
            dX = Eigen::Vector3d(d_point3d(gen),d_point3d(gen),d_point3d(gen));
        }
        else
        {
            dX.setZero();
        }
        map_points_noisy[id] = MapPoint(id,map_point.second.X + dX);

        for(const auto & ob : obs)
        {
            Eigen::Vector2d dp;
            if(cfg_noise.perturb_point2d)
            {
                dp = Eigen::Vector2d(d_point2d(gen),d_point2d(gen));
            }
            else
            {
                dp.setZero();
            }
            Eigen::Vector2d p_noisy = ob.second + dp;
            map_points_noisy[id].addObservation(ob.first,p_noisy);
        }
    }
    std::cout<<"add noise for points finished."<<std::endl;
}

void SimulationGenerator::Save(Simulation_Data& sfm_data, const std::string& outPath)
{
    if(access(outPath.c_str(),00) == -1)
    {
//        mkdir(outPath.c_str(),S_IRWXU);
        std::string cmd = "mkdir -p " + outPath;
        system(cmd.c_str());
    }
    std::string kp_outPath = outPath + "/keypoints/";
    if(access(kp_outPath.c_str(),00) == -1)
    {
//        mkdir(outPath.c_str(),S_IRWXU);
        std::string cmd = "mkdir -p " + kp_outPath;
        system(cmd.c_str());
    }

    MapPoints& map_points = sfm_data.map_points;
    KeyFrames& key_frames = sfm_data.key_frames;

    // save trajectory
    using namespace std;
    ofstream f(outPath + "/trajectory.txt");
    f<<fixed;
    f<<setprecision(7);
    double timestamp = 0.0;
    double deltaT = mpPoseGenerator->getDeltaT();
    for(auto & key_frame : key_frames)
    {
        Pose& p = key_frame.second.pose;

        Eigen::Quaterniond q = p.q.inverse();
        Eigen::Vector3d t = - (q * p.t);
        f<<timestamp<<" "<<t[0]<<" "<<t[1]<<" "<<t[2]<<" "
         <<q.coeffs()[0]<<" "<<q.coeffs()[1]<<" "<<q.coeffs()[2]<<" "<<q.coeffs()[3]<<endl;

        timestamp += deltaT;
    }
    f.close();
    std::cout<<"trajectory saved."<<std::endl;

    // save cameras
    std::vector<double> params;
    mpCamera->getParams(params);
    f.open(outPath + "/cameras.txt");
    f<<"# Camera list with one line of data per camera:"<<endl;
    f<<"#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]"<<endl;
    f<<"# Number of cameras: 1"<<endl;
    f<<"0 "<<mpCamera->getModelName()<<" "<<mpCamera->getWidth()<<" "<<mpCamera->getHeight();
    for(int i=0;i<params.size();i++)
    {
        f<<" "<<params[i];
    }
    f<<endl;
    f.close();
    cout<<"cameras saved."<<endl;

    // save images
    int count = 0;
    f.open(outPath + "/images.txt");
    f<<"# Image list with two lines of data per image:"<<endl;
    f<<"#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME"<<endl;
    f<<"#   POINTS2D[] as (X, Y, POINT3D_ID)"<<endl;
    f<<"# Number of images: "<<key_frames.size()<<", mean observations per image: "<<count/key_frames.size()<<endl;

    f<<fixed;
    f<<setprecision(7);
    int kf_count = 0;
    for(auto & key_frame : key_frames)
    {
        string img_name = "VIRTUAL_IMG_"+std::to_string(kf_count);
        Eigen::Quaterniond& q = key_frame.second.pose.q;
        Eigen::Vector3d& t = key_frame.second.pose.t;
        f<<key_frame.second.Id<<" "<<q.coeffs()[3]<<" "<<q.coeffs()[0]<<" "<<q.coeffs()[1]<<" "<<q.coeffs()[2]<<" "
         <<t[0]<<" "<<t[1]<<" "<<t[2]<<" "<<("0 " + img_name + ".JPG")<<endl;

        // saving KeyPoints
        ofstream kp_f(kp_outPath + img_name + ".feat");

        std::vector<unsigned int>& obs = key_frame.second.obs;
        bool flag = false;
        for(unsigned int point3d_id : obs)
        {
            Eigen::Vector2d& xp = map_points[point3d_id].obs[key_frame.second.Id];
            if(flag)
            {
                f<<" "<<xp[0]<<" "<<xp[1]<<" "<<point3d_id;
                kp_f<<endl<<xp[0]<<" "<<xp[1]<<" "<<1.0<<" "<<0.0;
            }
            else
            {
                f<<xp[0]<<" "<<xp[1]<<" "<<point3d_id;
                kp_f<<xp[0]<<" "<<xp[1]<<" "<<1.0<<" "<<0.0;
                flag = true;
            }
        }
        f<<endl;
        kp_f.close();
        kf_count++;
    }
    f.close();
    cout<<"images saved."<<endl;

    // save point3D
    f.open(outPath + "/points3D.txt");
    f<<"# 3D point list with one line of data per point:"<<endl;
    f<<"#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)"<<endl;
    f<<"# Number of points: "<<map_points.size()<<", mean track length: 3.3334"<<endl;

    f<<fixed;
    f<<setprecision(7);
    for(auto & map_point : map_points)
    {
        Eigen::Vector3d& X = map_point.second.X;
        f<<map_point.second.Id<<" "<<X[0]<<" "<<X[1]<<" "<<X[2]<<" 0 0 0 0.0";
        for(auto & kf : map_point.second.obs)
        {
            unsigned int image_id = kf.first;
            KeyFrame& key_frame = key_frames[image_id];
            int point2d_idx = -1;
            int idx = 0;
            for(auto & point_idx : key_frame.obs)
            {
                if(point_idx == map_point.second.Id)
                {
                    point2d_idx = idx;
                    break;
                }
                idx++;
            }
            f<<" "<<image_id<<" "<<point2d_idx;
        }
        f<<endl;
    }
    f.close();
    cout<<"points saved."<<endl;
}

}  // namespace generator