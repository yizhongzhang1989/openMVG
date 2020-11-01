//
// Created by v-xinli1 on 11/1/2020.
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

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <io.h>
#include <openMVG/sfm/sfm.hpp>
#include <nonFree/sift/SIFT_describer.hpp>
#include <dependencies/cereal/include/cereal/archives/json.hpp>


#include "random_cir.hpp"

class CsvReader_Simu
{
public:
    CsvReader_Simu() = delete;
    CsvReader_Simu( const std::string& filename, const char split = ' ' )
    {
        std::cout << "open file" << filename << std::endl;
        split_= split;
        csvInput_.open( filename );
        if( !csvInput_ ) throw std::runtime_error( "Error csv file dict: " + filename );

        std::string header;
        getline( csvInput_, header );

        data_ = std::vector<double>(8, 0);
        last_data_ = std::vector<double>(8, 0);
    }
    ~CsvReader_Simu()
    {
        if(csvInput_) csvInput_.close();
    }

    bool readline()
    {
        if( !csvInput_ ) throw std::runtime_error(" Not Set File to Read ");
        if( csvInput_.eof() ) return false;

        std::string line;
        getline(csvInput_, line);

        std::istringstream readStr( line );
        std::string part;

        for( int i= 0;i<8;++i  )
        {
            getline( readStr, part, split_ );
            data_[i] = std::strtod( part.c_str(), nullptr );
        }

        data_[0] *= 1e3;
        data_[0] = static_cast<long long int>(data_[0]);

        if( last_data_[0] != 0 )
        {
            if( (data_[0] - last_data_[0]) != 5 )
            {
                std::cout << line << std::endl;
                std::cout << data_[0] << " - " <<last_data_[0] << " != 5" << std::endl;
                return false;
            }
        }
        last_data_ = data_;
        return true;
    }

    std::vector<double> data_;
private:
    std::vector<double> last_data_;
    std::ifstream csvInput_;
    char split_;
};

struct KeyFrame
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    unsigned int Id_;
    Eigen::Matrix3d Rwc_;
    Eigen::Vector3d twc_;
    std::vector<unsigned long long int> obs_;
    KeyFrame()
            :Id_(-1)
    {
        obs_.clear();
    }
    KeyFrame(unsigned int id, const Eigen::Matrix3d& Rwc, const Eigen::Vector3d& twc)
            : Id_(id)
    {
        Rwc_ = Rwc;
        twc_ = twc;
        obs_.clear();
    }
    void addObservation(unsigned int PTId)
    {
        obs_.push_back(PTId);
    }
};

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace openMVG::features;
using namespace openMVG::matching;

std::vector<KeyFrame> KeyFrames;


int main(int argc, char **argv)
{
    CmdLine cmd;

    std::string simu_data_file;
    std::string sSfM_Data_dir;

    cmd.add( make_option('i', simu_data_file, "input_file") );
    cmd.add( make_option('o', sSfM_Data_dir, "output_dir") );

    cmd.process(argc, argv);

    CsvReader_Simu reader( simu_data_file );


    Eigen::Matrix3d Ric;
    Eigen::Vector3d tic;
    Ric << 0.0148655429818, -0.999880929698, 0.00414029679422,
            0.999557249008, 0.0149672133247, 0.025715529948,
            -0.0257744366974, 0.00375618835797, 0.999660727178;
    tic << -0.0216401454975, -0.064676986768, 0.00981073058949;

    Eigen::Quaterniond Qic(Ric);

    std::vector<Eigen::Vector3d> twcs;
    std::vector<Eigen::Quaterniond> Qwcs;
    while(reader.readline())
    {
        Eigen::Vector3d twi( reader.data_[1],reader.data_[2],reader.data_[3] );
        Eigen::Quaterniond Qwi( reader.data_[7],reader.data_[4],reader.data_[5],reader.data_[6] );
        Eigen::Vector3d twc = twi + Qwi*tic;
        Eigen::Quaterniond Qwc = Qwi * Qic;
        twcs.push_back(twc);
        Qwcs.push_back(Qwc);
    }

    Eigen::Vector3d center(0.,0.,0.);
    for(auto&twc:twcs)
    {
        center += twc;
    }
    center /= twcs.size();

    double max = 0., mean = 0.;

    for(auto&twc:twcs)
    {
        Eigen::Vector3d dis = twc - center;
        double distacne = dis.norm();
        max = max>distacne?max:distacne;
        mean += distacne;
    }
    mean /= twcs.size();

    std::cout << "max = " << max << std::endl;
    std::cout << "mean = " << mean << std::endl;


    int key_frame_index = 0;
    for( int i=0;i<twcs.size();++i )
    {
        if(i%20 == 0)
        {
            Eigen::Matrix3d Rwc = Qwcs[i].toRotationMatrix();
            KeyFrame kf( key_frame_index, Rwc, twcs[i] );
            KeyFrames.push_back(kf);
            key_frame_index++;
        }
    }

    std::vector<Eigen::Vector3d> Points;
    int z1 = 20;
    int z2=112;
    cRandom sita(z1,0.0);
    cRandom pesi(z2,0.0);

    const int num_point = 1e3;
    const double scale = 5.0;

    for(int i=0;i!=num_point;++i)
    {
        sita=my_random(pesi.seed);
        pesi=my_random(sita.seed);

        double u=2*sita.random-1.0;
        double v=2*pesi.random-1.0;

        double r2=pow(u,2)+pow(v,2);
        if(r2<1)
        {
            double x=2*u*sqrt(1-r2);
            double y=2*v*sqrt(1-r2);
            double z=1-2*r2;

            Eigen::Vector3d point(x, y, z);
            point *= scale;
            point += center;
            Points.push_back(point);
        }
    }


    std::vector< std::map<unsigned int, Eigen::Vector3d> > points_obs2ds;
    const double focal = 500;
    const double cx = 510;
    const double cy = 440;
    const double width = 1020;
    const double height = 880;
    int index_point = 0;
    for( auto&point:Points )
    {
        std::map<unsigned int, Eigen::Vector3d> obs2ds_;
        for( auto&kf:KeyFrames )
        {
            Eigen::Vector3d twc = kf.twc_;
            Eigen::Matrix3d Rwc = kf.Rwc_;
            Eigen::Matrix3d Rcw = Rwc.transpose();
            Eigen::Vector3d tcw = -Rcw*twc;
            Eigen::Vector3d point_cam = Rcw*point + tcw;

            if( point_cam(2) < 0 ) continue;
            point_cam /= point_cam(2);

            double pixle_x = focal*point_cam(0) + cx;
            double pixle_y = focal*point_cam(1) + cy;
            if( pixle_x < 0 || pixle_y < 0 ) continue;
            if( pixle_x >= width || pixle_y >= height ) continue;

            kf.addObservation(index_point);
            Eigen::Vector3d obs2d(pixle_x, pixle_y, 1);
            obs2ds_.insert( std::make_pair( kf.Id_, obs2d ) );
        }
        std::cout << "obs2ds_.size() = " << obs2ds_.size() << std::endl;
        points_obs2ds.push_back(obs2ds_);
        index_point++;
    }


    SfM_Data sfm_data;
    sfm_data.s_root_path = sSfM_Data_dir; // Setup main image root_path
    Views & views = sfm_data.views;
    Intrinsics & intrinsics = sfm_data.intrinsics;

    std::shared_ptr<IntrinsicBase> intrinsic;
    intrinsic = std::make_shared<Pinhole_Intrinsic>
            (width, height, focal, cx, cy);

    for(int i=0;i<KeyFrames.size();i++)
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
         stlplus::create_filespec( sSfM_Data_dir, "sfm_data.json" ).c_str(),
         ESfM_Data(VIEWS|INTRINSICS));

    std::cout << std::endl
              << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
              << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;


//    std::unique_ptr<Image_describer> image_describer;
//    const std::string sImage_describer = stlplus::create_filespec(sSfM_Data_dir, "image_describer", "json");
//    image_describer.reset(new SIFT_Image_describer
//                                  (SIFT_Image_describer::Params()));
//    image_describer->Set_configuration_preset(features::EDESCRIBER_PRESET(-1));
//    // Export the used Image_describer and region type for:
//    // - dynamic future regions computation and/or loading
//    {
//        std::ofstream stream(sImage_describer.c_str());
//        if (!stream.is_open())
//            return -1;
//
//        cereal::JSONOutputArchive archive(stream);
//        archive(cereal::make_nvp("image_describer", image_describer));
//        auto regionsType = image_describer->Allocate();
//        archive(cereal::make_nvp("regions_type", regionsType));
//    }


    // save key points
    int kf_count = 0;
    for(auto & key_frame : KeyFrames)
    {
        std::string img_name = "VIRTUAL_IMG_"+std::to_string(kf_count);
        std::ofstream kp_f(sSfM_Data_dir + "/" + img_name + ".feat");

        std::vector<uint64_t>& obs = key_frame.obs_;
        bool flag = false;
        for(unsigned int point3d_id : obs)
        {
            Eigen::Vector3d xp = points_obs2ds[point3d_id][key_frame.Id_];// map_points[point3d_id].obs[key_frame.second.Id];
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
    std::cout << "SavePoint suc" << std::endl;

    // save matches
    PairWiseMatches map_Matches;
    Pair_Set pairs;
    const size_t N = KeyFrames.size();
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

        auto kf1 = KeyFrames[I];

        for (unsigned int J : indexToCompare)
        {
            auto kf2 = KeyFrames[J];

            IndMatches vec_putatives_matches;
            {
                std::vector<uint64_t>& obs1 = kf1.obs_;
                std::vector<uint64_t>& obs2 = kf2.obs_;
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
                map_Matches.insert( { {I,J}, vec_putatives_matches } );
            }
        }
    }

    std::cout << "Generate match suc" << std::endl;

    if (!Save(map_Matches, std::string(sSfM_Data_dir + "/matches.f.bin")))
    {
        std::cerr
                << "Cannot save computed matches in: "
                << std::string(sSfM_Data_dir + "/matches.putative.txt");
    }

    return 0;
}