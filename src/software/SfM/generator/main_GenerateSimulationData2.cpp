#include <iostream>
#include <thread>
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
#include "PoseGeneratorSampling.h"
#include "SimulationGenerator.h"
#include "SurfaceSampler.h"
#include "TrajectorySampler.h"
#include "Utils.h"

#ifdef USE_PANGOLIN
#include "Visualizer.h"
typedef slam_visualization::Visualizer<Eigen::Vector3d, generator::InversePose, Eigen::aligned_allocator<Eigen::Vector3d>, Eigen::aligned_allocator<generator::InversePose>> GeneratorVisualizer;
#endif

using namespace generator;
using namespace std;

void SaveSfMData(std::string sImageDir, generator::Simulation_Data& simulationData, generator::CameraPinhole* pCam);
void SaveToImages(generator::Simulation_Data& sfm_data, const std::string& outPath, generator::CameraPinhole* pCam);
bool ParseExtrinsics(const std::string& str, std::vector<double>& params);

#ifdef USE_PANGOLIN
void DrawGrid(STLVector<Eigen::Vector3d>& up_vertices, STLVector<Eigen::Vector3d>& down_vertices)
{
    assert(up_vertices.size() == down_vertices.size());
    int N = up_vertices.size();

    glLineWidth(2);
    glBegin(GL_LINES);

    for(int i = 0; i < N - 1; i++)
    {
        Eigen::Vector3d& O = down_vertices[i];
        Eigen::Vector3d& A = down_vertices[i+1];
        Eigen::Vector3d& B = up_vertices[i+1];
        Eigen::Vector3d& C = up_vertices[i];

        glColor3f(1.0,0,1.0);
        glVertex3f(O[0],O[1],O[2]); glVertex3f(A[0],A[1],A[2]);
        glVertex3f(C[0],C[1],C[2]); glVertex3f(B[0],B[1],B[2]);

        glColor3f(0,1.0,1.0);
        glVertex3f(O[0],O[1],O[2]); glVertex3f(C[0],C[1],C[2]);
        glVertex3f(A[0],A[1],A[2]); glVertex3f(B[0],B[1],B[2]);
    }

    glEnd();
}
#endif

int main(int argc, char** argv)
{
    CmdLine cmd;

    std::string sTrajectoryObjFile;
    std::string sModelObjFile;
    std::string sOutDir;
    std::string sExtrinsics;
    double total_duration = 10.0;
    double f_cam = 30.0;
    int f_imu = 200;

    // generator config
    SimulationGenerator::SimulationConfig cfg;
    cfg.n_points = 1500;
    cfg.n_poses = -1;
    cfg.image_width = 640;
    cfg.image_height = 480;

    // required
    cmd.add(make_option('p', sTrajectoryObjFile, "trajectory_file"));
    cmd.add(make_option('m', sModelObjFile, "model_file"));
    cmd.add(make_option('o', sOutDir, "outdir"));
    // Optional
    cmd.add(make_option('t',total_duration,"duration"));
    cmd.add(make_option('e',sExtrinsics,"extrinsics"));
    cmd.add(make_option('n',cfg.n_points,"n_points"));
    cmd.add(make_option('f',f_cam,"f_cam"));
    cmd.add(make_option('F',f_imu,"f_imu"));

    try
    {
        if (argc == 1) throw std::string("Invalid command line parameter.");
        cmd.process(argc, argv);
    }
    catch (const std::string& s)
    {
        std::cerr << "Usage: " << argv[0] << '\n'
                  << "[-p|--trajectory_file] Obj file for generating trajectory \n"
                  << "[-m|--model_file] Obj file for generating structure (3D model) \n"
                  << "[-o|--outdir] output directory \n"
                  << "\n[Optional]\n"
                  << "[-t|--duration] total duration of trajectory in seconds, default is 10.0s \n"
                  << "[-e|--extrinsics] extrinsics from IMU to camera (T_cam_imu), format: (tx ty tz qx qy qz qw) \n"
                  << "[-n|--n_points] number of points in simulation \n"
                  << "[-f|--f_cam] frequency of camera (fps) \n"
                  << "[-F|--f_imu] frequency of imu (Hz), must be divisible for 1000 \n"
                  << std::endl;
        std::cerr << s << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << " You called : " <<std::endl
              << argv[0] << std::endl
              << "--trajectory_file " << sTrajectoryObjFile << std::endl
              << "--model_file " << sModelObjFile << std::endl
              << "--outdir "<< sOutDir << std::endl
              << "--duration " << total_duration << std::endl
              << "--extrinsics " << sExtrinsics << std::endl
              << "--n_points " << cfg.n_points << std::endl
              << "--f_cam " << f_cam << std::endl
              << "--f_imu " << f_imu << std::endl
              << std::endl;

    if(sOutDir.empty())
    {
        std::cerr<<"Output directory must be specified."<<std::endl;
        return EXIT_FAILURE;
    }
    Utils::check_and_create_dir(sOutDir);

    if(sTrajectoryObjFile.empty() || !stlplus::file_exists(sTrajectoryObjFile))
    {
        std::cerr<<"Invalid trajectory file: "<<sTrajectoryObjFile<<std::endl;
        return EXIT_FAILURE;
    }

    if(sModelObjFile.empty() || !stlplus::file_exists(sModelObjFile))
    {
        std::cerr<<"Invalid model file: "<<sModelObjFile<<std::endl;
        return EXIT_FAILURE;
    }

    if(1000 % f_imu)
    {
        std::cerr<<"frequency of IMU must be divisible for 1000."<<std::endl;
        return EXIT_FAILURE;
    }

    double T_cam = 1.0 / f_cam;
    int T_IMU = 1000 / f_imu;

    PoseGeneratorSampling g_pose(sTrajectoryObjFile,total_duration,T_cam,T_IMU,true);
    PointGenerator g_point(sModelObjFile);
    CameraPinhole cam(320,320,320,240,cfg.image_width,cfg.image_height);

    SimulationGenerator g_sim(&g_pose,&g_point,&cam);

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
    g_sim.SaveIMU(g_sim.getIMUMeasurements<PoseGeneratorSampling>(),openMVGPath + "/imu_data.csv");

#ifdef USE_PANGOLIN
    STLVector<Eigen::Vector3d> vertices;
    STLVector<InversePose> inv_poses = TrajectorySampler::SampleTrajectory(sTrajectoryObjFile,&vertices);
//    STLVector<InversePose> inv_poses_interpolated = TrajectorySampler::SampleTrajectory(sTrajectoryObjFile,10.0,0.05);
    const STLVector<InversePose>& inv_poses_interpolated = g_pose.getInversePoses();
    const STLVector<InversePose>& inv_poses_imu = g_pose.getInversePosesIMU();
    double deltaT_IMU = 5 * 1e-3;
    std::cout<<"Initial velocity: "<<((inv_poses_imu[1].p-inv_poses_imu[0].p)/deltaT_IMU).transpose()<<std::endl;

    int N = vertices.size() / 2;
    STLVector<Eigen::Vector3d> up_vertices, down_vertices;
    for(int i = 0; i < N; i++)
    {
        down_vertices.push_back(vertices[2 * i]);
        up_vertices.push_back(vertices[2 * i + 1]);
    }

    STLVector<Eigen::Vector3d> points;
    const MapPoints& map_points = sfm_data.map_points;
    for(const auto & map_point : map_points)
    {
        points.push_back(map_point.second.X);
    }

    // setup visualizer
    auto GetPoseCenter = [](const InversePose& pose)->Eigen::Vector3d
    {
        return pose.p;
    };
    auto Pose2Matrix = [](const InversePose& pose)->pangolin::OpenGlMatrix
    {
        Eigen::Matrix4d Twc;
        Twc << pose.q.toRotationMatrix(), pose.p, Eigen::Vector3d::Zero().transpose(), 1;
        return {Twc};
    };
    GeneratorVisualizer::VisualizerConfig vis_cfg;
    vis_cfg.GetPoseCenter = GetPoseCenter;
    vis_cfg.Pose2Matrix = Pose2Matrix;
    GeneratorVisualizer visualizer("Visualization",vis_cfg);

    pangolin::OpenGlRenderState* s_cam = visualizer.GetRenderState();
    pangolin::View* d_cam = visualizer.GetView();
    GeneratorVisualizer::ColorOptions color_opt;

    // create menu
    pangolin::CreatePanel("menu").SetBounds(0.0,1.0,0.0,pangolin::Attach::Pix(175));
    pangolin::Var<bool> interpolate_pose("menu.interpolate pose",false,true);

    while(!pangolin::ShouldQuit())
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        d_cam->Activate(*s_cam);
        glClearColor(1.0f,1.0f,1.0f,1.0f);

        visualizer.DrawPoints(points,color_opt.point_color);

        if(interpolate_pose)
        {
            visualizer.DrawPoses(inv_poses_interpolated,false,color_opt.pose_color);
            visualizer.DrawGraph(inv_poses_interpolated,color_opt.graph_color);
        }
        else
        {
            visualizer.DrawPoses(inv_poses,false,color_opt.pose_color);
            visualizer.DrawGraph(inv_poses,color_opt.graph_color);
        }

        DrawGrid(up_vertices,down_vertices);

        pangolin::FinishFrame();
    }

#endif

    return 0;
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

    if (!Save(map_Matches, std::string(sImageDir + "/matches.putative.bin")))
    {
        std::cerr
                << "Cannot save computed matches in: "
                << std::string(sImageDir + "/matches.putative.bin");
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