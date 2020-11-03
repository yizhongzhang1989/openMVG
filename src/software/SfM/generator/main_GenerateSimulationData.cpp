#include <iostream>

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

#include "PointGenerator.h"
#include "CameraPinhole.h"
#include "PoseGeneratorCircleSine.h"
#include "PoseGeneratorConstAcc.h"
#include "SimulationGenerator.h"
#include "SurfaceSampler.h"

void SaveSfMData(std::string sImageDir, generator::Simulation_Data& simulationData, generator::CameraPinhole* pCam);
void SaveToImages(generator::Simulation_Data& sfm_data, const std::string& outPath, generator::CameraPinhole* pCam);

#ifdef _WIN32
static std::string unix_path_to_win(const std::string& unix_format)
{
    std::string win_format(unix_format);
    for (int i = 0; i < win_format.length(); i++)
    {
        if (win_format[i] == '/')
        {
            win_format[i] = '\\';
        }
    }

    return win_format;
}
#endif

static void check_and_create_dir(const std::string& path)
{
#ifdef _WIN32
    if (_access(path.c_str(), 00) == -1)
    {
        std::string cmd = "mkdir " + unix_path_to_win(path);
        system(cmd.c_str());
    }
#elif __linux__
    if (access(path.c_str(), 00) == -1)
    {
        std::string cmd = "mkdir -p " + path;
        system(cmd.c_str());
    }
#endif
}

int main()
{
    using namespace generator;

    PointGenerator g_point(-25,25,-25,25,-25,25);
    PoseGeneratorCircleSine g_pose(5.0,0.1,1.0,0.1,0.005,true);
//    PoseGeneratorConstAcc g_pose(0.2,0.0,0.0,0.05,true,generator::PoseGeneratorConstAcc::FORWARD);
    CameraPinhole cam(320,320,320,240,640,480);

    // generation
    SimulationGenerator g_sim(&g_pose,&g_point,&cam);
    SimulationGenerator::SimulationConfig cfg;
    cfg.n_points = 1500;
    cfg.n_poses = 200;
    cfg.image_width = 640;
    cfg.image_height = 480;
    Simulation_Data sfm_data;
    g_sim.Generate(sfm_data,cfg);
    g_sim.Save(sfm_data,"groundtruth");
    SaveToImages(sfm_data,"images",&cam);

    // add noise
    Simulation_Data sfm_data_noisy;
    SimulationGenerator::NoiseConfig cfg_noise;
    cfg_noise.perturb_point2d = true;
    cfg_noise.perturb_point3d = true;
    cfg_noise.perturb_translation = true;
    cfg_noise.perturb_rotation = true;
    cfg_noise.stddev_point2d = 1.0;
    cfg_noise.stddev_point3d = 0.1;
    cfg_noise.stddev_translation = 0.1;
    cfg_noise.stddev_rotation = 0.01;
    g_sim.AddNoise(sfm_data,sfm_data_noisy,cfg_noise);
    g_sim.Save(sfm_data_noisy,"noisyData");

    SaveSfMData("openMVGFormat",sfm_data,&cam);

    g_sim.SaveIMU(g_sim.getIMUMeasurements<PoseGeneratorCircleSine>(),"openMVGFormat/imu_data.csv");
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

void SaveToImages(generator::Simulation_Data& sfm_data, const std::string& imagePath, generator::CameraPinhole* pCam)
{
    using namespace openMVG::image;

    check_and_create_dir(imagePath);

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
        std::string img_name = imagePath + "/VIRTUAL_IMG_"+std::to_string(kf_count) + ".JPG";
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