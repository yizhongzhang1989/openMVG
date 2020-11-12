#include <iostream>
#include <memory>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/features/akaze/image_describer_akaze_io.hpp"
#include "nonFree/sift/SIFT_describer_io.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <cereal/archives/json.hpp>
#include <cereal/details/helpers.hpp>
#include "openMVG/features/regions_factory_io.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"
#include "openMVG/image/image_io.hpp"
#include "third_party/cmdLine/cmdLine.h"

#include "PointGenerator.h"
#include "CameraPinhole.h"
#include "PoseGeneratorCircleSine.h"
#include "PoseGeneratorConstAcc.h"
#include "PoseGeneratorLine.h"
#include "SimulationGenerator.h"
#include "Utils.h"

#ifdef USE_PANGOLIN
#include "Visualizer.h"
typedef slam_visualization::Visualizer<Eigen::Vector3d, generator::InversePose, Eigen::aligned_allocator<Eigen::Vector3d>, Eigen::aligned_allocator<generator::InversePose>> GeneratorVisualizer;
#endif

void SaveSfMData(std::string sImageDir, generator::Simulation_Data& simulationData, generator::CameraPinhole* pCam);
void SaveToImages(generator::Simulation_Data& sfm_data, const std::string& outPath, generator::CameraPinhole* pCam);
bool ParseExtrinsics(const std::string& str, std::vector<double>& params);

int main(int argc, char** argv)
{
    using namespace generator;

    CmdLine cmd;

    std::string sObjFile = "";
    std::string sOutDir = "";
    std::string sExtrinsics = "";

    // generator config
    SimulationGenerator::SimulationConfig cfg;
    cfg.n_points = 1500;
    cfg.n_poses = 200;
    cfg.image_width = 640;
    cfg.image_height = 480;

    // bound
    Bound bound;
    bound.max_x = 25.0, bound.min_x = -25.0;
    bound.max_y = 25.0, bound.min_y = -25.0;
    bound.max_z = 25.0, bound.min_z = -25.0;

    enum TrajectoryType
    {
        CIRCLE_SINE = 0,
        CONST_ACC = 1,
        LINE = 2
    };
    TrajectoryType trajectory_type = CIRCLE_SINE;
    int trajectory_type_aux = 0;

    // required
    cmd.add( make_option('o', sOutDir, "outdir") );
    // Optional
    cmd.add(make_option('i', sObjFile, "obj_file"));
    cmd.add(make_option('n',cfg.n_points,"n_points"));
    cmd.add(make_option('N',cfg.n_poses,"n_poses"));
    cmd.add(make_option('x',bound.min_x,"min_x"));
    cmd.add(make_option('X',bound.max_x,"max_x"));
    cmd.add(make_option('y',bound.min_y,"min_y"));
    cmd.add(make_option('Y',bound.max_y,"max_y"));
    cmd.add(make_option('z',bound.min_z,"min_z"));
    cmd.add(make_option('Z',bound.max_z,"max_z"));
    cmd.add(make_option('p',trajectory_type_aux,"trajectory_type"));
    cmd.add(make_option('e',sExtrinsics,"extrinsics"));

    try
    {
        if (argc == 1) throw std::string("Invalid command line parameter.");
        cmd.process(argc, argv);
    }
    catch (const std::string& s)
    {
        std::cerr << "Usage: " << argv[0] << '\n'
                  << "[-o|--outdir path] \n"
                  << "\n[Optional]\n"
                  << "[-i|--obj_file] OBJ file representing 3D model (for point sampling) \n"
                  << "[-n|--n_points] number of points in simulation \n"
                  << "[-N|--n_poses] number of poses in simulation \n"
                  << "[-x|--min_x] minimum x value (range of points) \n"
                  << "[-X|--max_x] maximum x value (range of points) \n"
                  << "[-y|--min_y] minimum y value (range of points) \n"
                  << "[-Y|--max_y] maximum y value (range of points) \n"
                  << "[-z|--min_z] minimum z value (range of points) \n"
                  << "[-Z|--max_z] maximum z value (range of points) \n"
                  << "[-p|--trajectory_type] trajectory type, currently support : [0|CIRCLE_SINE, 1|CONST_ACC, 2|lINE] \n"
                  << "[-e|--extrinsics] extrinsics from IMU to camera (T_cam_imu), format: (tx ty tz qx qy qz qw) \n"
                  << std::endl;
        std::cerr << s << std::endl;
        return EXIT_FAILURE;
    }

    trajectory_type = static_cast<TrajectoryType>(trajectory_type_aux);

    std::cout << " You called : " <<std::endl
              << argv[0] << std::endl
              << "--obj_file " << sObjFile << std::endl
              << "--outdir " << sOutDir << std::endl
              << "--n_points " << cfg.n_points << std::endl
              << "--n_poses " << cfg.n_poses << std::endl
              << "--min_x " << bound.min_x << std::endl
              << "--max_x " << bound.max_x << std::endl
              << "--min_y " << bound.min_y << std::endl
              << "--max_y " << bound.max_y << std::endl
              << "--min_z " << bound.min_z << std::endl
              << "--max_z " << bound.max_z << std::endl
              << "--trajectory_type " << trajectory_type << std::endl
              << "--extrinsics " << sExtrinsics <<std::endl
              << std::endl;

    if(sOutDir.empty())
    {
        std::cerr<<"Output directory must be specified."<<std::endl;
        return EXIT_FAILURE;
    }
    Utils::check_and_create_dir(sOutDir);

    // create pose generator
    std::shared_ptr<PoseGeneratorBase<Pose, Eigen::aligned_allocator<Pose>>> g_pose;
    switch(trajectory_type)
    {
        case CIRCLE_SINE:
        {
            g_pose = std::make_shared<PoseGeneratorCircleSine>(5.0,0.2,0.1,1.0,50,5,true);
        }break;
        case CONST_ACC:
        {
            g_pose = std::make_shared<PoseGeneratorConstAcc>(0.2,0.0,0.0,50,5,true,generator::PoseGeneratorConstAcc::FORWARD);
        }break;
        case LINE:
        {
            g_pose = std::make_shared<PoseGeneratorLine>(50,5,true);
        }break;
        default:
        {
            std::cerr<<"Unrecognized trajectory type : "<<trajectory_type<<std::endl;
            return EXIT_FAILURE;
        }
    }

    // create point generator
    std::shared_ptr<PointGeneratorBase<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>>> g_point;
    if(sObjFile.empty())
    {
        g_point = std::make_shared<PointGenerator>(bound.min_x,bound.max_x,bound.min_y,bound.max_y,bound.min_z,bound.max_z);
    }
    else
    {
        g_point = std::make_shared<PointGenerator>(sObjFile);
    }

    // create camera
    CameraPinhole cam(320,320,320,240,cfg.image_width,cfg.image_height);

    SimulationGenerator g_sim(g_pose.get(),g_point.get(),&cam);

    if(!sExtrinsics.empty())
    {
        std::vector<double> params;
        if(ParseExtrinsics(sExtrinsics,params))
        {
            Pose T_cam_imu;
            T_cam_imu.t = Eigen::Vector3d(params[0],params[1],params[2]);
            T_cam_imu.q = Eigen::Quaterniond(params[6],params[3],params[4],params[5]);
            T_cam_imu.q.normalize();
            g_sim.setExtrinsics(T_cam_imu);
        }
        else
        {
            std::cerr<<"failed to parse extrinsics: "<<sExtrinsics<<std::endl;
            return EXIT_FAILURE;
        }
    }

    Simulation_Data sfm_data;
    g_sim.Generate(sfm_data,cfg);

    // save as colmap format
    std::string groundtruthPath = sOutDir + "/groundtruth/";
    Utils::check_and_create_dir(groundtruthPath);
    g_sim.Save(sfm_data,groundtruthPath);
    // save to images (for visualization)
    std::string imagePath = sOutDir + "/images/";
    Utils::check_and_create_dir(imagePath);
    SaveToImages(sfm_data,imagePath,&cam);
    // save as openMVG format
    std::string openMVGPath = sOutDir + "/openMVGFormat/";
    Utils::check_and_create_dir(openMVGPath);
    SaveSfMData(openMVGPath,sfm_data,&cam);
    // save IMU
    switch(trajectory_type)
    {
        case CIRCLE_SINE:
        {
            g_sim.SaveIMU(g_sim.getIMUMeasurements<PoseGeneratorCircleSine>(),openMVGPath + "/imu_data.csv");
        }break;
        case CONST_ACC:
        {
            g_sim.SaveIMU(g_sim.getIMUMeasurements<PoseGeneratorConstAcc>(),openMVGPath + "/imu_data.csv");
        }break;
        case LINE:
        {
            g_sim.SaveIMU(g_sim.getIMUMeasurements<PoseGeneratorLine>(),openMVGPath + "/imu_data.csv");
        }break;
        default:
        {
            std::cerr<<"Unrecognized trajectory type : "<<trajectory_type<<std::endl;
            return EXIT_FAILURE;
        }
    }

//    // add noise
//    Simulation_Data sfm_data_noisy;
//    SimulationGenerator::NoiseConfig cfg_noise;
//    cfg_noise.perturb_point2d = true;
//    cfg_noise.perturb_point3d = true;
//    cfg_noise.perturb_translation = true;
//    cfg_noise.perturb_rotation = true;
//    cfg_noise.stddev_point2d = 1.0;
//    cfg_noise.stddev_point3d = 0.1;
//    cfg_noise.stddev_translation = 0.1;
//    cfg_noise.stddev_rotation = 0.01;
//    g_sim.AddNoise(sfm_data,sfm_data_noisy,cfg_noise);
//    g_sim.Save(sfm_data_noisy,"noisyData");
}

void SaveSfMData(std::string sImageDir, generator::Simulation_Data& simulationData, generator::CameraPinhole* pCam)
{
    using namespace openMVG;
    using namespace openMVG::sfm;
    using namespace openMVG::cameras;
    using namespace openMVG::features;
    using namespace openMVG::matching;

    if(access(sImageDir.c_str(),00) == -1)
    {
//        mkdir(outPath.c_str(),S_IRWXU);
        std::string cmd = "mkdir -p " + sImageDir;
        system(cmd.c_str());
    }

    // Configure an empty scene with Views and their corresponding cameras
    SfM_Data sfm_data;
    sfm_data.s_root_path = sImageDir; // Setup main image root_path
    Views & views = sfm_data.views;
    Intrinsics & intrinsics = sfm_data.intrinsics;

    // Expected properties for each image
    double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;
    std::vector<double> params;
    pCam->getParams(params);
    width = pCam->getWidth();
    height = pCam->getHeight();
    focal = params[0];
    ppx = params[2];
    ppy = params[3];

    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic;
    intrinsic = std::make_shared<Pinhole_Intrinsic>
            (width, height, focal, ppx, ppy);

    for(int i=0;i<simulationData.key_frames.size();i++)
    {
        View v("VIRTUAL_IMG_"+std::to_string(views.size())+".JPG", views.size(), views.size(), views.size(), width, height);

        // Add intrinsic related to the image (if any)
        if (intrinsic == nullptr)
        {
            //Since the view have invalid intrinsic data
            // (export the view, with an invalid intrinsic field value)
            v.id_intrinsic = UndefinedIndexT;
        }
        else
        {
            // Add the defined intrinsic to the sfm_container
            intrinsics[v.id_intrinsic] = intrinsic;
        }
        // Add the view to the sfm_container
        views[v.id_view] = std::make_shared<View>(v);
    }
    GroupSharedIntrinsics(sfm_data);

    // Store SfM_Data views & intrinsic data
    Save(sfm_data,
        stlplus::create_filespec( sImageDir, "sfm_data.json" ).c_str(),
        ESfM_Data(VIEWS|INTRINSICS));

    std::cout << std::endl
              << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
              << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;

    std::unique_ptr<Image_describer> image_describer;
    const std::string sImage_describer = stlplus::create_filespec(sImageDir, "image_describer", "json");
    image_describer.reset(new SIFT_Image_describer
                                  (SIFT_Image_describer::Params()));
    image_describer->Set_configuration_preset(features::EDESCRIBER_PRESET(-1));
    // Export the used Image_describer and region type for:
    // - dynamic future regions computation and/or loading
    {
        std::ofstream stream(sImage_describer.c_str());
        if (!stream.is_open())
            return;

        cereal::JSONOutputArchive archive(stream);
        archive(cereal::make_nvp("image_describer", image_describer));
        auto regionsType = image_describer->Allocate();
        archive(cereal::make_nvp("regions_type", regionsType));
    }

    generator::MapPoints& map_points = simulationData.map_points;
    generator::KeyFrames& key_frames = simulationData.key_frames;

    // save key points
    int kf_count = 0;
    for(auto & key_frame : key_frames)
    {
        std::string img_name = "VIRTUAL_IMG_"+std::to_string(kf_count);
        std::ofstream kp_f(sImageDir + "/" + img_name + ".feat");

        std::vector<unsigned int>& obs = key_frame.second.obs;
        bool flag = false;
        for(unsigned int point3d_id : obs)
        {
            Eigen::Vector2d& xp = map_points[point3d_id].obs[key_frame.second.Id];
            if(flag)
            {
                kp_f<<std::endl<<xp[0]<<" "<<xp[1]<<" "<<1.0<<" "<<0.0;
            }
            else
            {
                kp_f<<xp[0]<<" "<<xp[1]<<" "<<1.0<<" "<<0.0;
                flag = true;
            }
        }
        kp_f.close();

        kf_count++;
    }

    // save matches
    PairWiseMatches map_Matches;
    Pair_Set pairs;
    const size_t N = key_frames.size();
    for(int i=0;i<N;i++)
    {
        for(int j=i+1;j<N;j++)
        {
            pairs.insert({i,j});
        }
    }
    // Sort pairs according the first index to minimize the MatcherT build operations
    using Map_vectorT = std::map<IndexT, std::vector<IndexT>>;
    Map_vectorT map_Pairs;
    for (const auto & pair_it : pairs)
    {
        map_Pairs[pair_it.first].push_back(pair_it.second);
    }
    // Perform matching between all the pairs
    for (const auto & pairs_it : map_Pairs)
    {
        const IndexT I = pairs_it.first;
        const auto & indexToCompare = pairs_it.second;

        generator::KeyFrame& kf1 = key_frames[I];

        for (int j = 0; j < static_cast<int>(indexToCompare.size()); ++j)
        {
            const IndexT J = indexToCompare[j];

            generator::KeyFrame& kf2 = key_frames[J];

            IndMatches vec_putatives_matches;
            {
                std::vector<unsigned int>& obs1 = kf1.obs;
                std::vector<unsigned int>& obs2 = kf2.obs;
                for(int il = 0; il < obs1.size(); il++)
                {
                    for(int ir = 0; ir < obs2.size(); ir++)
                    {
                        if(obs1[il] == obs2[ir])
                        {
                            vec_putatives_matches.emplace_back(il,ir);
                            break;
                        }
                    }
                }
            }
            if (!vec_putatives_matches.empty())
            {
                map_Matches.insert( { {I,J}, std::move(vec_putatives_matches) } );
            }
        }
    }

    if (!Save(map_Matches, std::string(sImageDir + "/matches.putative.txt")))
    {
        std::cerr
                << "Cannot save computed matches in: "
                << std::string(sImageDir + "/matches.putative.txt");
    }
}

void SaveToImages(generator::Simulation_Data& sfm_data, const std::string& outPath, generator::CameraPinhole* pCam)
{
    using namespace openMVG::image;

    Utils::check_and_create_dir(outPath);

    int width = pCam->getWidth();
    int height = pCam->getHeight();
    generator::MapPoints& map_points = sfm_data.map_points;
    generator::KeyFrames& key_frames = sfm_data.key_frames;

    auto assignColor = [](Image<RGBColor>& img, RGBColor& color, int u, int v, int w)
    {
        int startU = u - w;
        int endU = u + w;
        int startV = v - w;
        int endV = v + w;

        if(startU < 0)
            startU = 0;
        if(endU >= img.Width())
            endU = img.Width()-1;
        if(startV < 0)
            startV = 0;
        if(endV >= img.Height())
            endV = img.Height()-1;

        for(int i = startU; i < endU; i++)
        {
            for(int j = startV; j < endV; j++)
            {
                img(j,i) = color;
            }
        }
    };

    int kf_count = 0;
    for(auto & key_frame : key_frames)
    {
        std::string img_name = outPath + "/VIRTUAL_IMG_"+std::to_string(kf_count) + ".JPG";
        Image<RGBColor> img(width,height,true,RGBColor(255,255,255));
        std::vector<unsigned int>& obs = key_frame.second.obs;
        for(unsigned int point3d_id : obs)
        {
            Eigen::Vector2d& xp = map_points[point3d_id].obs[key_frame.second.Id];
            generator::Color& color = map_points[point3d_id].color;

            int u = int(xp.x());
            int v = int(xp.y());
//            img(v,u) = RGBColor(color.r,color.g,color.b);
            RGBColor rgb(color.r,color.g,color.b);
            assignColor(img,rgb,u,v,2);
        }
        WriteImage(img_name.c_str(),img);
        kf_count++;
    }
}

bool ParseExtrinsics(const std::string& str, std::vector<double>& params)
{
    std::stringstream ss;
    ss<<str;
    double d;
    params.clear();
    for(int i = 0; i < 7 && !ss.eof(); i++)
    {
        ss>>d;
        params.push_back(d);
    }
    return params.size() == 7;
}