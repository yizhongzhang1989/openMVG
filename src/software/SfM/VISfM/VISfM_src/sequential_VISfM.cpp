//
// Created by root on 10/13/20.
//

#include "sequential_VISfM.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#include <ceres/types.h>
#include <functional>
#include <iostream>
#include <utility>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

    using namespace openMVG::cameras;
    using namespace openMVG::geometry;
    using namespace openMVG::matching;

    SequentialVISfMReconstructionEngine::SequentialVISfMReconstructionEngine(
            const SfM_Data & sfm_data,
            const std::string &soutDirectory,
            const std::string &sloggingFile)
            : ReconstructionEngine(sfm_data, soutDirectory),
              sLogging_file_(sloggingFile),
              initial_pair_(0,0),
              cam_type_(EINTRINSIC(PINHOLE_CAMERA_RADIAL3))
    {
        if (!sLogging_file_.empty())
        {
         // setup HTML logger
         html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("SequentialReconstructionEngine SFM report.");
         html_doc_stream_->pushInfo(
                 htmlDocument::htmlMarkup("h1", std::string("SequentialSfMReconstructionEngine")));
         html_doc_stream_->pushInfo("<hr>");

         html_doc_stream_->pushInfo( "Dataset info:");
         html_doc_stream_->pushInfo( "Views count: " +
                                     htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
        }
        // Init remaining image list
        for (Views::const_iterator itV = sfm_data.GetViews().begin();
          itV != sfm_data.GetViews().end(); ++itV)
        {
            set_remaining_view_id_.insert(itV->second->id_view);
        }
    }

    SequentialVISfMReconstructionEngine::~SequentialVISfMReconstructionEngine()
    {
        if (!sLogging_file_.empty())
        {
            // Save the reconstruction Log
            std::ofstream htmlFileStream(sLogging_file_.c_str());
            htmlFileStream << html_doc_stream_->getDoc();
        }
    }

    void SequentialVISfMReconstructionEngine::SetFeaturesProvider(Features_Provider * provider)
    {
        features_provider_ = provider;
    }

    void SequentialVISfMReconstructionEngine::SetMatchesProvider(Matches_Provider * provider)
    {
        matches_provider_ = provider;
    }

    void SequentialVISfMReconstructionEngine::SetTimeStamp(std::vector<IndexT> &times)
    {
        assert( times.size() == sfm_data_.views.size() );
        int index = 0;
        for( auto&Id_View:sfm_data_.views )
        {
            assert( index < times.size() );
            auto Id = Id_View.second->id_pose;
            sfm_data_.timestamps.emplace( Id, times[index++] );
        }
    }

    void SequentialVISfMReconstructionEngine::SetIMUDataset(std::shared_ptr<IMU_Dataset> imudataset_)
    {
        sfm_data_.imu_dataset = std::move(imudataset_);
    }

    bool SequentialVISfMReconstructionEngine::Process()
    {
        // Compute robust Resection of remaining images
        // - group of images will be selected and resection + scene completion will be tried
        size_t resectionGroupIndex = 0;
        std::vector<uint32_t> vec_possible_resection_indexes;
        while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
        {
            bool bImageAdded = false;
            // Add images to the 3D reconstruction
            for (const auto & iter : vec_possible_resection_indexes)
            {
                bImageAdded |= Resection(iter);
                set_remaining_view_id_.erase(iter);
            }

            if (bImageAdded)
            {
                // Scene logging as ply for visual debug
                std::ostringstream os;
                os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
                Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

                // Perform BA until all point are under the given precision
                do
                {
                    update_state_speed();
                    update_imu_time();
                    update_imu_inte();
                    BundleAdjustment_optimizi_only_IMU();
                    std::cout << "end BundleAdjustment_optimizi_only_IMU" << std::endl;
                    BundleAdjustmentWithIMU();
                    std::cout << "end BundleAdjustmentWithIMU" << std::endl;
                }
                while (badTrackRejector(4.0, 50));
                eraseUnstablePosesAndObservations(sfm_data_);
            }
            ++resectionGroupIndex;
        }
        // Ensure there is no remaining outliers
        if (badTrackRejector(4.0, 0))
        {
            eraseUnstablePosesAndObservations(sfm_data_);
        }

        //-- Reconstruction done.
        //-- Display some statistics
        std::cout << "\n\n-------------------------------" << "\n"
                  << "-- Structure from Motion (statistics):\n"
                  << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
                  << " from " << sfm_data_.GetViews().size() << " input images.\n"
                  << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
                  << "-------------------------------" << "\n";

        Histogram<double> h;
        ComputeResidualsHistogram(&h);
        std::cout << "\nHistogram of residuals:\n" << h.ToString() << std::endl;

        if (!sLogging_file_.empty())
        {
            using namespace htmlDocument;
            std::ostringstream os;
            os << "Structure from Motion process finished.";
            html_doc_stream_->pushInfo("<hr>");
            html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

            os.str("");
            os << "-------------------------------" << "<br>"
               << "-- Structure from Motion (statistics):<br>"
               << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
               << " from " <<sfm_data_.GetViews().size() << " input images.<br>"
               << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "<br>"
               << "-------------------------------" << "<br>";
            html_doc_stream_->pushInfo(os.str());

            html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

            const std::vector<double> xBin = h.GetXbinsValue();
            const auto range = autoJSXGraphViewport<double>(xBin, h.GetHist());

            htmlDocument::JSXGraphWrapper jsxGraph;
            jsxGraph.init("3DtoImageResiduals",600,300);
            jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
            jsxGraph.UnsuspendUpdate();
            jsxGraph.setViewport(range);
            jsxGraph.close();
            html_doc_stream_->pushInfo(jsxGraph.toStr());
        }

        {
            {
                auto it_pose =  sfm_data_.poses.begin();
                auto tw0 = it_pose->second.center();
                auto it0 = it_pose->first;
                while( std::next(it_pose) != sfm_data_.poses.end() ) it_pose++;
                auto tw1 = it_pose->second.center();
                auto it1 = it_pose->first;
                std::cout << "it0 = " << it0 << std::endl;
                std::cout << "it1 = " << it1 << std::endl;
                std::cout << "t01 = " << (tw1 - tw0).norm() << std::endl;
            }
        }
        return true;
    }

    bool SequentialVISfMReconstructionEngine::Process_onlyvisual()
    {
        return false;
    }

    bool SequentialVISfMReconstructionEngine::VI_Init()
    {
        //-------------------
        //-- Incremental reconstruction
        //-------------------

        if (!InitLandmarkTracks())
            return false;

        // Initial pair choice
        if (initial_pair_ == Pair(0,0))
        {
            // Initial pair must be choice already

            // TODO xinli select check
            if (!AutomaticInitialPairChoice(initial_pair_))
            {
                std::cerr << "Cannot find a valid initial pair" << std::endl;
                // Cannot find a valid initial pair, try to set it by hand?
                if (!ChooseInitialPair(initial_pair_))
                {
                    return false;
                }
            }

            std::cout << "---------------------------------------\n"
                      << "initial_pair_.first = " << initial_pair_.first << "\n"
                      <<  "initial_pair_.second = " << initial_pair_.second << "\n"
                      << "---------------------------------------" << std::endl;

            // get start index and end index of vi init
            {
                const int len_size_2 = 15;
                IndexT len_2 = (initial_pair_.second - initial_pair_.first) / 2;
                IndexT inital_len_2 = len_size_2 > len_2 ? len_size_2 : len_2;
                IndexT center = initial_pair_.first + len_2;
                IndexT left, right;
                if( center < inital_len_2 )
                    left = 0;
                else
                    left = center - inital_len_2;
                right = center + inital_len_2;
//                if( left < 0 )
//                {
////                    right -= left;
//                    left = 0;
////                    assert( right >= sfm_data_.GetViews().size() );
//                }
                if( right >= sfm_data_.GetViews().size())
                {
//                    int over = right - sfm_data_.GetViews().size() + 1;
//                    left -= over;
                    right = sfm_data_.GetViews().size() -1;
//                    assert( left < 0 );
                }

                // Init remaining image list
                for (Views::const_iterator itV = sfm_data_.GetViews().begin();
                     itV != sfm_data_.GetViews().end(); ++itV)
                {
                    if( itV->second->id_view >= left && itV->second->id_view <= right )
                    {
                        set_remaining_view_id_vi_init_.insert(itV->second->id_view);
                        set_remaining_view_id_.erase(itV->second->id_view);
//                        set_remaining_view_id_.erase(v_id);
                    }
                }
            }
        }
        // Else a starting pair was already initialized before

        // Initial pair Essential Matrix and [R|t] estimation.
        if (!MakeInitialPair3D(initial_pair_))
        {
            std::cerr << "initial pair solve error" << std::endl;
            return false;
        }

        // Compute robust Resection of remaining images
        // - group of images will be selected and resection + scene completion will be tried
        size_t resectionGroupIndex = 0;
        std::vector<uint32_t> vec_possible_resection_indexes;
        while (FindImagesWithPossibleResection_VIinit(vec_possible_resection_indexes))
        {
            bool bImageAdded = false;
            // Add images to the 3D reconstruction
            for (const auto & iter : vec_possible_resection_indexes)
            {
                bImageAdded |= Resection(iter);
                set_remaining_view_id_vi_init_.erase(iter);
            }

            if (bImageAdded)
            {
                // Scene logging as ply for visual debug
                std::ostringstream os;
                os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
                Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

                // Perform BA until all point are under the given precision
                do
                {
                    BundleAdjustmentVisualInit();
//                    BundleAdjustment();
                }
                while (badTrackRejector(4.0, 50));
                eraseUnstablePosesAndObservations(sfm_data_);
            }
            ++resectionGroupIndex;
        }
        // Ensure there is no remaining outliers
        if (badTrackRejector(4.0, 0))
        {
            eraseUnstablePosesAndObservations(sfm_data_);
        }

        //-- Reconstruction done.
        //-- Display some statistics
        std::cout << "\n\n-------------------------------" << "\n"
                  << "-- Structure from Motion (statistics):\n"
                  << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
                  << " from " << sfm_data_.GetViews().size() << " input images.\n"
                  << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
                  << "-------------------------------" << "\n";

        Histogram<double> h;
        ComputeResidualsHistogram(&h);
        std::cout << "\nHistogram of residuals:\n" << h.ToString() << std::endl;

        if (!sLogging_file_.empty())
        {
            using namespace htmlDocument;
            std::ostringstream os;
            os << "Structure from Motion process finished.";
            html_doc_stream_->pushInfo("<hr>");
            html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

            os.str("");
            os << "-------------------------------" << "<br>"
               << "-- Structure from Motion (statistics):<br>"
               << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
               << " from " <<sfm_data_.GetViews().size() << " input images.<br>"
               << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "<br>"
               << "-------------------------------" << "<br>";
            html_doc_stream_->pushInfo(os.str());

            html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

            const std::vector<double> xBin = h.GetXbinsValue();
            const auto range = autoJSXGraphViewport<double>(xBin, h.GetHist());

            htmlDocument::JSXGraphWrapper jsxGraph;
            jsxGraph.init("3DtoImageResiduals",600,300);
            jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
            jsxGraph.UnsuspendUpdate();
            jsxGraph.setViewport(range);
            jsxGraph.close();
            html_doc_stream_->pushInfo(jsxGraph.toStr());
        }


        for(  auto& v_id:set_remaining_view_id_vi_init_ )
        {
            set_remaining_view_id_.insert(v_id);
        }
        return true;
    }

    void SequentialVISfMReconstructionEngine::rota_pose()
    {
        for( auto& it_pose:sfm_data_.poses )
        {
            Mat3 Rcw = it_pose.second.rotation();
            Rcw = sfm_data_.IG_Ric * Rcw;
            it_pose.second.SetRoation(Rcw);
        }
//        for( auto&landmark:sfm_data_.structure )
//        {
//            landmark.second.X = ( Ric.transpose() * landmark.second.X );
//        }
    }

    void SequentialVISfMReconstructionEngine::recover_g_s(const Eigen::Vector3d& correct_g, const Eigen::VectorXd& speeds_scale)
    {
        Eigen::Matrix3d R0 = Utility::g2R(correct_g);

        Mat3 Rw0 = sfm_data_.poses.begin()->second.rotation().transpose();

        double yaw = Utility::R2ypr(R0 * Rw0).x();
        R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
//        std::cout << "correct_g = " << correct_g.transpose() << std::endl;
//        correct_g = R0 * correct_g;
//        std::cout << "correct_g = " << correct_g.transpose() << std::endl;
        //Matrix3d rot_diff = R0 * Rs[0].transpose();
        Eigen::Quaterniond rot_diff( R0 );
        VIstaticParm::G_ = R0 * VIstaticParm::G_;

//        for( auto& it_pose:sfm_data_.poses )
//        {
//            Mat3 Rcw = it_pose.second.rotation();
//            Rcw = sfm_data_.IG_Ric.transpose() * Rcw;
//            it_pose.second.SetRoation(Rcw);
//        }

        /*{
            for (auto & structure_landmark_it : sfm_data_.structure)
            {
                const Observations & obs = structure_landmark_it.second.obs;

                Vec3 Point3D =structure_landmark_it.second.X;

                for (const auto & obs_it : obs)
                {
                    // Build the residual block corresponding to the track observation:
                    const View * view = sfm_data_.views.at(obs_it.first).get();
                    Vec3 tiw = sfm_data_.poses[view->id_pose].translation();
                    Mat3 Riw = sfm_data_.poses[view->id_pose].rotation();

                    Mat3 Rci = sfm_data_.IG_Ric.transpose();
                    Vec3 tci = -Rci * sfm_data_.IG_tic;

                    std::vector<double> intri = sfm_data_.intrinsics[0]->getParams();
                    Vec3 PointI = Riw * Point3D + tiw;
                    Vec3 PointC = Rci * PointI;
                    PointC /= PointC(2);

                    const double &focal = intri[0];
                    const double &principal_point_x = intri[1];
                    const double &principal_point_y = intri[2];
                    const double &k1 = intri[3];
                    const double &k2 = intri[4];
                    const double &k3 = intri[5];

                    // Apply distortion (xd,yd) = disto(x_u,y_u)
                    const double r2 = PointC.head(2).squaredNorm();
                    const double r4 = r2 * r2;
                    const double r6 = r4 * r2;
                    const double r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);

                    Vec2 residuals;
                    residuals << principal_point_x + (PointC.x() * r_coeff) * focal - obs_it.second.x.x(),
                            principal_point_y + (PointC.y() * r_coeff) * focal - obs_it.second.x.y();

                    std::cout << "residuals = " << residuals.transpose() << std::endl;

                }
            }
        }*/

        std::cout <<"-=------=-=-=-=-=-=-=-=" << std::endl;
        double correct_scale = (speeds_scale.tail<1>())(0) / 100.;


        int kv = -1;
        auto pose_i = sfm_data_.poses.begin();
        for (; pose_i != sfm_data_.poses.end(); pose_i++)
        {
            kv++;
            auto poseId = pose_i->first;
            auto Rwi = pose_i->second.rotation().transpose();
            Vec3 speedV3d = Rwi * speeds_scale.segment<3>(kv * 3);
            IMU_Speed speed(speedV3d);
            sfm_data_.Speeds.insert( std::make_pair(poseId, speed) );
        }

        for( auto& it_pose:sfm_data_.poses )
        {
            Mat3 Riw = it_pose.second.rotation();
            Riw = Riw * rot_diff.inverse();
            it_pose.second.SetRoation(Riw);

            Vec3 twi = it_pose.second.center();
            twi = correct_scale*(rot_diff * twi) - Riw.transpose() * sfm_data_.IG_tic;
            it_pose.second.SetCenter(twi);

            Mat3 Rcw = sfm_data_.IG_Ric.transpose() * Riw;
            Vec3 twc = twi + Riw.transpose() * sfm_data_.IG_tic;
            it_pose.second.SetRoation(Rcw);
            it_pose.second.SetCenter(twc);
        }

        for( auto&landmark:sfm_data_.structure )
        {
            landmark.second.X = correct_scale*(rot_diff * landmark.second.X);
        }

        for( auto&speed:sfm_data_.Speeds )
        {
            speed.second.speed_ = rot_diff * speed.second.speed_;
        }

        /*{
            for (auto & structure_landmark_it : sfm_data_.structure)
            {
                const Observations & obs = structure_landmark_it.second.obs;

                Vec3 Point3D =structure_landmark_it.second.X;

                for (const auto & obs_it : obs)
                {
                    // Build the residual block corresponding to the track observation:
                    const View * view = sfm_data_.views.at(obs_it.first).get();
                    Vec3 tiw = sfm_data_.poses[view->id_pose].translation();
                    Mat3 Riw = sfm_data_.poses[view->id_pose].rotation();

                    Mat3 Rci = sfm_data_.IG_Ric.transpose();
                    Vec3 tci = -Rci * sfm_data_.IG_tic;

                    std::vector<double> intri = sfm_data_.intrinsics[0]->getParams();
                    Vec3 PointI = Riw * Point3D + tiw;
                    Vec3 PointC = Rci * PointI + tci;
                    PointC /= PointC(2);

                    const double &focal = intri[0];
                    const double &principal_point_x = intri[1];
                    const double &principal_point_y = intri[2];
                    const double &k1 = intri[3];
                    const double &k2 = intri[4];
                    const double &k3 = intri[5];

                    // Apply distortion (xd,yd) = disto(x_u,y_u)
                    const double r2 = PointC.head(2).squaredNorm();
                    const double r4 = r2 * r2;
                    const double r6 = r4 * r2;
                    const double r_coeff = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);

                    Vec2 residuals;
                    residuals << principal_point_x + (PointC.x() * r_coeff) * focal - obs_it.second.x.x(),
                            principal_point_y + (PointC.y() * r_coeff) * focal - obs_it.second.x.y();

                    std::cout << "residuals = " << residuals.transpose() << std::endl;

                }
            }
        }*/
    }

    void SequentialVISfMReconstructionEngine::update_imu_inte()
    {
        for( auto & id_imubase:sfm_data_.imus )
        {
            IndexT t0 = id_imubase.second.t0_;
            IndexT t1 = id_imubase.second.t1_;
            std::vector< IndexT > times;
            std::vector<Vec3> accs;
            std::vector<Vec3> gyrs;

//            std::cout << "start GetMeasure" << std::endl;
            std::tie( times, accs, gyrs ) = sfm_data_.imu_dataset->GetMeasure(t0, t1);
//            std::cout << "start integrate" << std::endl;
            id_imubase.second.integrate(accs, gyrs, times);
        }
    }

    void SequentialVISfMReconstructionEngine::update_imu_time()
    {
//        sfm_data_.imus.clear();
        auto it_pose = sfm_data_.poses.begin();
        IndexT last_t = sfm_data_.timestamps.at( it_pose->first );
        it_pose++;
        while( it_pose != sfm_data_.poses.end() )
        {
            IndexT cur_t = sfm_data_.timestamps.at( it_pose->first );
            auto id_pose = it_pose->first;
            if( sfm_data_.imus.count( id_pose ) == 0 )
            {
                //todo xinli ba bg
//                std::cout << "start insert imu" << std::endl;

//                std::cout << "start construction imubase_ptr" << std::endl;
                IMU_InteBase imubase_ptr(last_t, cur_t);
//                std::shared_ptr<IMU_InteBase> imubase_ptr = std::make_shared<IMU_InteBase>(last_t, cur_t);
//                std::cout << "end construction imubase_ptr" << std::endl;
                sfm_data_.imus.insert( std::make_pair( id_pose,  imubase_ptr) );


//                std::cout << "start construction IMU_Speed" << std::endl;
//                IMU_Speed speed(Eigen::Vector3d(0.,0.,0.));
//                std::cout << "end construction IMU_Speed" << std::endl;
//                sfm_data_.Speeds.insert( std::make_pair(id_pose, speed) );
//                std::cout << "end insert imu" << std::endl;
            }
            else
            {
                sfm_data_.imus[id_pose].change_time(last_t, cur_t);
            }
            last_t = cur_t;
            it_pose++;
        }
    }

    void SequentialVISfMReconstructionEngine::update_state_speed()
    {
        auto it_pose = sfm_data_.poses.begin();
        while( it_pose != sfm_data_.poses.end() )
        {
            IndexT cur_t = sfm_data_.timestamps.at( it_pose->first );
            auto id_pose = it_pose->first;
            if( sfm_data_.Speeds.count( id_pose ) == 0 )
            {
//                std::cout << "start construction IMU_Speed" << std::endl;
                IMU_Speed speed(Eigen::Vector3d(0.,0.,0.)); // TODO xinli
//                std::cout << "end construction IMU_Speed" << std::endl;
                sfm_data_.Speeds.insert( std::make_pair(id_pose, speed) );
//                std::cout << "end insert imu" << std::endl;
            }
            it_pose++;
        }
    }

    bool SequentialVISfMReconstructionEngine::check_imu_observibility()
    {
        //check imu observibility
        {
            Eigen::Vector3d sum_g;
            double sum_i = 0.;
            for( const auto&imu:sfm_data_.imus )
            {
                auto dt = imu.second.sum_dt_;
                if( dt == 0 ) continue;
                Eigen::Vector3d tmp_g = imu.second.delta_v_ / dt;
                sum_g += tmp_g;
                sum_i += 1.;
            }

            Eigen::Vector3d aver_g;
            aver_g = sum_g * 1.0 / sum_i;
            double var = 0;
            for( const auto&imu:sfm_data_.imus )
            {
                auto dt = imu.second.sum_dt_;
                if( dt == 0 ) continue;
                Eigen::Vector3d tmp_g = imu.second.delta_v_ / dt;
                var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            }

            var = std::sqrt(var / sum_i);
            //ROS_WARN("IMU variation %f!", var);
            std::cout << "IMU variation " << var << std::endl;
            if(var < 0.25)
            {
                std::cout << "IMU excitation not enouth!" << std::endl;
                return false;
            }
            return true;
        }
    }

    bool SequentialVISfMReconstructionEngine::VI_align()
    {
        std::cout << "start VI align" << std::endl;
        std::cout << "start update_imu_time" << std::endl;
        update_imu_time();
        std::cout << "start update_imu_inte" << std::endl;
        update_imu_inte();
        std::cout << "start rota_pose" << std::endl;
        update_state_speed();


        {
            check_imu_observibility();
        }



        rota_pose();
        Eigen::VectorXd speeds_scale;
        Eigen::Vector3d correct_g;

        std::cout << "start solveGyroscopeBias" << std::endl;
        solveGyroscopeBias();
        std::cout << "end solveGyroscopeBias" << std::endl;

        std::cout << "start solve_vgs" << std::endl;
        if( !solve_vgs(speeds_scale, correct_g) ) return false;
        std::cout << "end solve_vgs" << std::endl;

        // init pre not very good
        {
            auto it_pose = sfm_data_.poses.begin();
            auto tw0 = it_pose->second.center();
            auto it0 = it_pose->first;
            while( std::next(it_pose) != sfm_data_.poses.end() ) it_pose++;
            auto tw1 = it_pose->second.center();
            auto it1 = it_pose->first;
            std::cout << "it0 = " << it0 << std::endl;
            std::cout << "it1 = " << it1 << std::endl;
            std::cout << "t01 = " << (tw1 - tw0).norm() << std::endl;
        }

        std::cout << "start recover_g_s" << std::endl;
        recover_g_s( correct_g, speeds_scale );
        std::cout << "end recover_g_s" << std::endl;

        {
            auto it_pose = sfm_data_.poses.begin();
            auto tw0 = it_pose->second.center();
            auto it0 = it_pose->first;
            while( std::next(it_pose) != sfm_data_.poses.end() ) it_pose++;
            auto tw1 = it_pose->second.center();
            auto it1 = it_pose->first;
            std::cout << "it0 = " << it0 << std::endl;
            std::cout << "it1 = " << it1 << std::endl;
            std::cout << "t01 = " << (tw1 - tw0).norm() << std::endl;
        }

//
//        BundleAdjustmentWithIMU();
//        BundleAdjustmentWithIMU();
        BundleAdjustment();
//        BundleAdjustment();
//        BundleAdjustmentWithIMU();
//        BundleAdjustmentWithIMU();
//        BundleAdjustment();
//        BundleAdjustment();

//        {
//            auto it_pose = sfm_data_.poses.begin();
//            auto tw0 = it_pose->second.center();
//            auto it0 = it_pose->first;
//            while( std::next(it_pose) != sfm_data_.poses.end() ) it_pose++;
//            auto tw1 = it_pose->second.center();
//            auto it1 = it_pose->first;
//            std::cout << "it0 = " << it0 << std::endl;
//            std::cout << "it1 = " << it1 << std::endl;
//            std::cout << "t01 = " << (tw1 - tw0).norm() << std::endl;
//        }

        return true;
    }

    void SequentialVISfMReconstructionEngine::solveGyroscopeBias()
    {
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        Eigen::Vector3d delta_bg;
        A.setZero();
        b.setZero();

        auto pose_it_i = sfm_data_.poses.begin();
        auto pose_it_j = sfm_data_.poses.begin(); pose_it_j++;
        for( ; pose_it_j != sfm_data_.poses.end(); pose_it_i++, pose_it_j++ )
        {
            Mat3 tmp_A;
            tmp_A.setZero();
            Eigen::Vector3d tmp_b;
            tmp_b.setZero();

            Mat3 Riw = pose_it_i->second.rotation();
//            Mat3 Rwi = Riw.transpose();
            Mat3 Rjw = pose_it_j->second.rotation();
            Mat3 Rwj = Rjw.transpose();

            auto imu_inte = sfm_data_.imus.at( pose_it_j->first );

            Eigen::Quaterniond q_ij(Riw * Rwj);
            tmp_A = imu_inte.GetBg();
            tmp_b = 2 * ( imu_inte.delta_q_.inverse() * q_ij).vec();
            A += tmp_A.transpose() * tmp_A;
            b += tmp_A.transpose() * tmp_b;


//            Eigen::Quaterniond delta = imu_inte.delta_q_.inverse() * q_ij;
//            std::cout << "delta = " << delta.coeffs().transpose() << std::endl;
        }
        delta_bg = A.ldlt().solve(b);

        for( auto&it_imu:sfm_data_.imus )
        {
            it_imu.second.repropagate( Eigen::Vector3d::Zero(), delta_bg );
        }
        std::cout << "delta_bg = " << delta_bg.transpose() << std::endl;
    }

    bool SequentialVISfMReconstructionEngine::solve_vgs(Eigen::VectorXd& speeds_scale, Eigen::Vector3d &correct_g)
    {
        int n_state = static_cast<int>(sfm_data_.imus.size()+1) * 3 + 3 + 1;
        Eigen::MatrixXd A{n_state, n_state};
        A.setZero();
        Eigen::VectorXd b{n_state};
        b.setZero();

        int index = 0;
        auto pose_it_i = sfm_data_.poses.begin();
        auto pose_it_j = sfm_data_.poses.begin(); pose_it_j++;
        for( ; pose_it_j != sfm_data_.poses.end(); pose_it_i++, pose_it_j++, index++ )
        {
            Eigen::MatrixXd tmp_A(6, 10);
            tmp_A.setZero();
            Eigen::VectorXd tmp_b(6);
            tmp_b.setZero();

            Mat3 Riw = pose_it_i->second.rotation();
            Mat3 Rwi = Riw.transpose();
            Mat3 Rjw = pose_it_j->second.rotation();
            Mat3 Rwj = Rjw.transpose();
            Vec3 twi = pose_it_i->second.center();
            Vec3 twj = pose_it_j->second.center();

            auto imu_inte = sfm_data_.imus.at( pose_it_j->first );

            double dt = imu_inte.sum_dt_;
            tmp_A.block<3, 3>(0, 0) = -dt * Eigen::Matrix3d::Identity();
            tmp_A.block<3, 3>(0, 6) = Rwi.transpose() * dt * dt / 2 * Eigen::Matrix3d::Identity();
            tmp_A.block<3, 1>(0, 9) = Rwi.transpose() * (twj - twi) / 100.;
            tmp_b.block<3, 1>(0, 0) = imu_inte.delta_p_ + Rwi.transpose() * Rwj * sfm_data_.IG_tic - sfm_data_.IG_tic;
            //cout << "delta_p   " << frame_j->second.pre_integratfion->delta_p.transpose() << endl;
            tmp_A.block<3, 3>(3, 0) = -Eigen::Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = Rwi.transpose() * Rwj;
            tmp_A.block<3, 3>(3, 6) = Rwi.transpose() * dt * Eigen::Matrix3d::Identity();
            tmp_b.block<3, 1>(3, 0) = imu_inte.delta_v_;
            //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

            Eigen::Matrix<double, 6, 6> cov_inv = Eigen::Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            Eigen::MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            Eigen::VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(index * 3, index * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(index * 3) += r_b.head<6>();

            A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
            b.tail<4>() += r_b.tail<4>();

            A.block<6, 4>(index * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
            A.block<4, 6>(n_state - 4, index * 3) += r_A.bottomLeftCorner<4, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        speeds_scale = A.ldlt().solve(b);
//        std::cout << "speeds_scale.size() = " << speeds_scale.size() << std::endl;
        double correct_scale = speeds_scale(n_state - 1) / 100.;
        correct_g = speeds_scale.segment<3>(n_state - 4);
        std::cout << "correct_scale = " << correct_scale << std::endl;
        std::cout << "correct_g nrom = " << correct_g.norm() << std::endl << "correct_g =  " << correct_g.transpose() << std::endl;

        if(std::fabs(correct_g.norm() - VIstaticParm::G_.norm()) > 2.0 || correct_scale < 0)
        {
            return false;
        }
        std::cout << "correct_scale = " << correct_scale << std::endl;
        std::cout << "correct_g nrom = " << correct_g.norm() << std::endl << "correct_g =  " << correct_g.transpose() << std::endl;

        RefineGravity(speeds_scale, correct_g);

        correct_scale = (speeds_scale.tail<1>())(0) / 100;
        std::cout << "correct_scale = " << correct_scale << std::endl;
        std::cout << "correct_g nrom = " << correct_g.norm() << std::endl << "correct_g =  " << correct_g.transpose() << std::endl;
        VIstaticParm::G_ = correct_g;
        return true;
    }

    Eigen::MatrixXd SequentialVISfMReconstructionEngine::TangentBasis(Eigen::Vector3d &g0)
    {
        Eigen::Vector3d b, c;
        Eigen::Vector3d a = g0.normalized();
        Eigen::Vector3d tmp(0, 0, 1);
        if(a == tmp)
            tmp << 1, 0, 0;
        b = (tmp - a * (a.transpose() * tmp)).normalized();
        c = a.cross(b);
        Eigen::MatrixXd bc(3, 2);
        bc.block<3, 1>(0, 0) = b;
        bc.block<3, 1>(0, 1) = c;
        return bc;
    }

    void SequentialVISfMReconstructionEngine::RefineGravity(Eigen::VectorXd& speeds_scale, Eigen::Vector3d &correct_g)
    {
        Eigen::Vector3d g0 = correct_g.normalized() * VIstaticParm::G_.norm();
        Eigen::Vector3d lx, ly;

        int n_state = static_cast<int>(sfm_data_.imus.size()+1) * 3 + 2 + 1;

        Eigen::MatrixXd A{n_state, n_state};
        A.setZero();
        Eigen::VectorXd b{n_state};
        b.setZero();
        for(int k = 0; k < 4; k++)
        {
            Eigen::MatrixXd lxly(3, 2);
            lxly = TangentBasis(g0);
            int index = 0;
            auto pose_it_i = sfm_data_.poses.begin();
            auto pose_it_j = sfm_data_.poses.begin(); pose_it_j++;
            for( ; pose_it_j != sfm_data_.poses.end(); pose_it_i++, pose_it_j++, index++ )
            {
                Eigen::MatrixXd tmp_A(6, 9);
                tmp_A.setZero();
                Eigen::VectorXd tmp_b(6);
                tmp_b.setZero();

                Mat3 Riw = pose_it_i->second.rotation();
                Mat3 Rwi = Riw.transpose();
                Mat3 Rjw = pose_it_j->second.rotation();
                Mat3 Rwj = Rjw.transpose();
                Vec3 twi = pose_it_i->second.center();
                Vec3 twj = pose_it_j->second.center();

                auto imu_inte = sfm_data_.imus.at( pose_it_j->first );
                double dt = imu_inte.sum_dt_;

                tmp_A.block<3, 3>(0, 0) = -dt * Eigen::Matrix3d::Identity();
                tmp_A.block<3, 2>(0, 6) = Riw * dt * dt / 2 * Eigen::Matrix3d::Identity() * lxly;
                tmp_A.block<3, 1>(0, 8) = Riw * (twj - twi) / 100.0;
                tmp_b.block<3, 1>(0, 0) = imu_inte.delta_p_ + Riw * Rwj * sfm_data_.IG_tic - sfm_data_.IG_tic - Riw * dt * dt / 2 * g0;

                tmp_A.block<3, 3>(3, 0) = - Eigen::Matrix3d::Identity();
                tmp_A.block<3, 3>(3, 3) = Riw * Rwj;
                tmp_A.block<3, 2>(3, 6) = Riw * dt * Eigen::Matrix3d::Identity() * lxly;
                tmp_b.block<3, 1>(3, 0) = imu_inte.delta_v_- Riw * dt * Eigen::Matrix3d::Identity() * g0;


                Eigen::Matrix<double, 6, 6> cov_inv = Eigen::Matrix<double, 6, 6>::Zero();
                //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
                //MatrixXd cov_inv = cov.inverse();
                cov_inv.setIdentity();

                Eigen::MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
                Eigen::VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

                A.block<6, 6>(index * 3, index * 3) += r_A.topLeftCorner<6, 6>();
                b.segment<6>(index * 3) += r_b.head<6>();

                A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
                b.tail<3>() += r_b.tail<3>();

                A.block<6, 3>(index * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
                A.block<3, 6>(n_state - 3, index * 3) += r_A.bottomLeftCorner<3, 6>();
            }
            A = A * 1000.0;
            b = b * 1000.0;
            speeds_scale = A.ldlt().solve(b);
            Eigen::VectorXd dg = speeds_scale.segment<2>(n_state - 3);
            g0 = (g0 + lxly * dg).normalized() * VIstaticParm::G_.norm();
            double correct_scale = speeds_scale(n_state-1) / 100.;
            std::cout << "refine correct_scale = " << correct_scale << std::endl;
            //double s = x(n_state - 1);
        }
        correct_g = g0;
    }

    bool SequentialVISfMReconstructionEngine::BundleAdjustmentVisualInit()
    {
        Bundle_Adjustment_Ceres::BA_Ceres_options options;
        if ( sfm_data_.GetPoses().size() > 100 &&
             (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
                )
            // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
        {
            options.preconditioner_type_ = ceres::JACOBI;
            options.linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
        else
        {
            options.linear_solver_type_ = ceres::DENSE_SCHUR;
        }
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
                ( ReconstructionEngine::intrinsic_refinement_options_,
                  Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
                  Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
                  Control_Point_Parameter(),
                  this->b_use_motion_prior_
                );
        return bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
    }

    bool SequentialVISfMReconstructionEngine::BundleAdjustment()
    {
        Bundle_Adjustment_IMU_Ceres::BA_Ceres_options options;
        if ( sfm_data_.GetPoses().size() > 100 &&
             (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
                )
            // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
        {
            options.preconditioner_type_ = ceres::JACOBI;
            options.linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
        else
        {
            options.linear_solver_type_ = ceres::DENSE_SCHUR;
        }
        Bundle_Adjustment_IMU_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
                ( ReconstructionEngine::intrinsic_refinement_options_,
                  Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
                  Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
                  Control_Point_Parameter(),
                  this->b_use_motion_prior_
                );
        return bundle_adjustment_obj.Adjust_onlyvisual(sfm_data_, ba_refine_options);
    }

    /// Bundle adjustment to refine Structure; Motion and Intrinsics
    bool SequentialVISfMReconstructionEngine::BundleAdjustmentWithIMU()
    {
        Bundle_Adjustment_IMU_Ceres::BA_Ceres_options options;
        if ( sfm_data_.GetPoses().size() > 100 &&
             (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
                )
            // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
        {
            options.preconditioner_type_ = ceres::JACOBI;
            options.linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
        else
        {
            options.linear_solver_type_ = ceres::DENSE_SCHUR;
        }
        Bundle_Adjustment_IMU_Ceres bundle_adjustment_obj(options);
//        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
                ( ReconstructionEngine::intrinsic_refinement_options_,
                  Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
                  Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
                  Control_Point_Parameter(),
                  this->b_use_motion_prior_
                );
        return bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
    }

    bool SequentialVISfMReconstructionEngine::BundleAdjustment_optimizi_only_IMU()
    {
        Bundle_Adjustment_IMU_Ceres::BA_Ceres_options options;
        if ( sfm_data_.GetPoses().size() > 100 &&
             (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
              ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
                )
            // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
        {
            options.preconditioner_type_ = ceres::JACOBI;
            options.linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
        else
        {
            options.linear_solver_type_ = ceres::DENSE_SCHUR;
        }
        Bundle_Adjustment_IMU_Ceres bundle_adjustment_obj(options);
//        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
                ( ReconstructionEngine::intrinsic_refinement_options_,
                  Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
                  Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
                  Control_Point_Parameter(),
                  this->b_use_motion_prior_
                );
        return bundle_adjustment_obj.Adjust_onlyIMU(sfm_data_, ba_refine_options);
    }

    /**
     * @brief Discard tracks with too large residual error
     *
     * Remove observation/tracks that have:
     *  - too large residual error
     *  - too small angular value
     *
     * @return True if more than 'count' outliers have been removed.
     */
    bool SequentialVISfMReconstructionEngine::badTrackRejector(double dPrecision, size_t count)
    {
        const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(sfm_data_, dPrecision, 2);
        const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(sfm_data_, 2.0);
        std::cout << "nbOutliers_residualErr = " << nbOutliers_residualErr << "\n"
        << "dPrecision = " << dPrecision << "\n"
        << "nbOutliers_angleErr = " << nbOutliers_angleErr << "\n"
        << "count = " << count << std::endl;

        return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
    }

    bool SequentialVISfMReconstructionEngine::InitLandmarkTracks()
    {
        // Compute tracks from matches
        tracks::TracksBuilder tracksBuilder;

        {
            // List of features matches for each couple of images
            const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
            std::cout << "\n" << "Track building" << std::endl;

            tracksBuilder.Build(map_Matches);
            std::cout << "\n" << "Track filtering" << std::endl;
            tracksBuilder.Filter();
            std::cout << "\n" << "Track export to internal struct" << std::endl;
            //-- Build tracks with STL compliant type :
            tracksBuilder.ExportToSTL(map_tracks_);

            std::cout << "\n" << "Track stats" << std::endl;
            {
                std::ostringstream osTrack;
                //-- Display stats :
                //    - number of images
                //    - number of tracks
                std::set<uint32_t> set_imagesId;
                tracks::TracksUtilsMap::ImageIdInTracks(map_tracks_, set_imagesId);
                osTrack << "------------------" << "\n"
                        << "-- Tracks Stats --" << "\n"
                        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
                        << " Images Id: " << "\n";
                std::copy(set_imagesId.begin(),
                          set_imagesId.end(),
                          std::ostream_iterator<uint32_t>(osTrack, ", "));
                osTrack << "\n------------------" << "\n";

                std::map<uint32_t, uint32_t> map_Occurence_TrackLength;
                tracks::TracksUtilsMap::TracksLength(map_tracks_, map_Occurence_TrackLength);
                osTrack << "TrackLength, Occurrence" << "\n";
                for (const auto & it : map_Occurence_TrackLength)  {
                    osTrack << "\t" << it.first << "\t" << it.second << "\n";
                }
                osTrack << "\n";
                std::cout << osTrack.str();
            }
        }
        // Initialize the shared track visibility helper
        shared_track_visibility_helper_.reset(new openMVG::tracks::SharedTrackVisibilityHelper(map_tracks_));
        return map_tracks_.size() > 0;
    }

    double SequentialVISfMReconstructionEngine::ComputeResidualsHistogram(Histogram<double> * histo)
    {
        // Collect residuals for each observation
        std::vector<float> vec_residuals;
        vec_residuals.reserve(sfm_data_.structure.size());
        for (const auto & landmark_entry : sfm_data_.GetLandmarks())
        {
            const Observations & obs = landmark_entry.second.obs;
            for (const auto & observation : obs)
            {
                const View * view = sfm_data_.GetViews().find(observation.first)->second.get();
                const Pose3 pose = sfm_data_.GetPoseOrDie(view);
                const auto intrinsic = sfm_data_.GetIntrinsics().find(view->id_intrinsic)->second;
                const Vec2 residual = intrinsic->residual(pose(landmark_entry.second.X), observation.second.x);
                vec_residuals.emplace_back( std::abs(residual(0)) );
                vec_residuals.emplace_back( std::abs(residual(1)) );
            }
        }
        // Display statistics
        if (vec_residuals.size() > 1)
        {
            float dMin, dMax, dMean, dMedian;
            minMaxMeanMedian<float>(vec_residuals.cbegin(), vec_residuals.cend(),
                                    dMin, dMax, dMean, dMedian);
            if (histo)  {
                *histo = Histogram<double>(dMin, dMax, 10);
                histo->Add(vec_residuals.cbegin(), vec_residuals.cend());
            }

            std::cout << std::endl << std::endl;
            std::cout << std::endl
                      << "SequentialSfMReconstructionEngine::ComputeResidualsMSE." << "\n"
                      << "\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size() << std::endl
                      << "\t-- Residual min:\t" << dMin << std::endl
                      << "\t-- Residual median:\t" << dMedian << std::endl
                      << "\t-- Residual max:\t "  << dMax << std::endl
                      << "\t-- Residual mean:\t " << dMean << std::endl;

            return dMean;
        }
        return -1.0;
    }

    /// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
    bool SequentialVISfMReconstructionEngine::MakeInitialPair3D(const Pair & current_pair)
    {
        // Compute robust Essential matrix for ImageId [I,J]
        // use min max to have I < J
        const uint32_t
                I = std::min(current_pair.first, current_pair.second),
                J = std::max(current_pair.first, current_pair.second);

        if (sfm_data_.GetViews().count(I) == 0 ||
                sfm_data_.GetViews().count(J) == 0)
        {
            return false;
        }
        // a. Assert we have valid cameras
        const View
                * view_I = sfm_data_.GetViews().at(I).get(),
                * view_J = sfm_data_.GetViews().at(J).get();
        const Intrinsics::const_iterator
                iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
                iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

        if (iterIntrinsic_I == sfm_data_.GetIntrinsics().end() ||
            iterIntrinsic_J == sfm_data_.GetIntrinsics().end() )
        {
            return false;
        }

        const auto
                * cam_I = iterIntrinsic_I->second.get(),
                * cam_J = iterIntrinsic_J->second.get();
        if (!cam_I || !cam_J)
        {
            return false;
        }

        // b. Get common features between the two view
        // use the track to have a more dense match correspondence set
        openMVG::tracks::STLMAPTracks map_tracksCommon;
        shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

        //-- Copy point to arrays
        const size_t n = map_tracksCommon.size();
        Mat xI(2,n), xJ(2,n);
        uint32_t cptIndex = 0;
        for (const auto & track_iter : map_tracksCommon)
        {
            auto iter = track_iter.second.cbegin();
            const uint32_t
                    i = iter->second,
                    j = (++iter)->second;

            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
            xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
            ++cptIndex;
        }

        // c. Robust estimation of the relative pose
        RelativePose_Info relativePose_info;

        const std::pair<size_t, size_t>
                imageSize_I(cam_I->w(), cam_I->h()),
                imageSize_J(cam_J->w(), cam_J->h());

        if (!robustRelativePose(
                cam_I, cam_J, xI, xJ, relativePose_info, imageSize_I, imageSize_J, 4096))
        {
            std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
                      << std::endl;
            return false;
        }
        std::cout << "A-Contrario initial pair residual: "
                  << relativePose_info.found_residual_precision << std::endl;
        // Bound min precision at 1 pix.
        relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

        const bool bRefine_using_BA = true;
        if (bRefine_using_BA)
        {
            // Refine the defined scene
            SfM_Data tiny_scene;
            tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
            tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
            tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
            tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

            // Init poses
            const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
            const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

            // Init structure
            Landmarks & landmarks = tiny_scene.structure;

            for (const auto & track_iterator : map_tracksCommon)
            {
                // Get corresponding points
                auto iter = track_iterator.second.cbegin();
                const uint32_t
                        i = iter->second,
                        j = (++iter)->second;

                const Vec2
                        x1 = features_provider_->feats_per_view[I][i].coords().cast<double>(),
                        x2 = features_provider_->feats_per_view[J][j].coords().cast<double>();

                Vec3 X;
                if (Triangulate2View(
                        Pose_I.rotation(),
                        Pose_I.translation(),
                        (*cam_I)(cam_I->get_ud_pixel(x1)),
                        Pose_J.rotation(),
                        Pose_J.translation(),
                        (*cam_J)(cam_J->get_ud_pixel(x2)),
                        X,
                        triangulation_method_))
                {
                    Observations obs;
                    obs[view_I->id_view] = Observation(x1, i);
                    obs[view_J->id_view] = Observation(x2, j);
                    landmarks[track_iterator.first].obs = std::move(obs);
                    landmarks[track_iterator.first].X = X;
                }
            }
            Save(tiny_scene, stlplus::create_filespec(sOut_directory_, "initialPair.ply"), ESfM_Data(ALL));

            // - refine only Structure and Rotations & translations (keep intrinsic constant)
            Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
            options.linear_solver_type_ = ceres::DENSE_SCHUR;
            Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
            if (!bundle_adjustment_obj.Adjust(tiny_scene,
                                              Optimize_Options
                                                      (
                                                              Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
                                                              Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
                                                              Structure_Parameter_Type::ADJUST_ALL) // Adjust structure
            )
                    )
            {
                return false;
            }

            // Save computed data
            const Pose3 pose_I = sfm_data_.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
            const Pose3 pose_J = sfm_data_.poses[view_J->id_pose] = tiny_scene.poses[view_J->id_pose];
            map_ACThreshold_.insert({I, relativePose_info.found_residual_precision});
            map_ACThreshold_.insert({J, relativePose_info.found_residual_precision});
            set_remaining_view_id_.erase(view_I->id_view);
            set_remaining_view_id_.erase(view_J->id_view);

            // List inliers and save them
            for (const auto & landmark_entry : tiny_scene.GetLandmarks())
            {
                const IndexT trackId = landmark_entry.first;
                const Landmark & landmark = landmark_entry.second;
                const Observations & obs = landmark.obs;
                Observations::const_iterator
                        iterObs_xI = obs.find(view_I->id_view),
                        iterObs_xJ = obs.find(view_J->id_view);

                const Observation & ob_xI = iterObs_xI->second;
                const Observation & ob_xJ = iterObs_xJ->second;
                const Vec2
                        ob_xI_ud = cam_I->get_ud_pixel(ob_xI.x),
                        ob_xJ_ud = cam_J->get_ud_pixel(ob_xJ.x);

                const double angle = AngleBetweenRay(
                        pose_I, cam_I, pose_J, cam_J, ob_xI_ud, ob_xJ_ud);
                const Vec2 residual_I = cam_I->residual(pose_I(landmark.X), ob_xI.x);
                const Vec2 residual_J = cam_J->residual(pose_J(landmark.X), ob_xJ.x);
                if (angle > 2.0 &&
                    CheiralityTest((*cam_I)(ob_xI_ud), pose_I,
                                   (*cam_J)(ob_xJ_ud), pose_J,
                                   landmark.X) &&
                    residual_I.norm() < relativePose_info.found_residual_precision &&
                    residual_J.norm() < relativePose_info.found_residual_precision)
                {
                    sfm_data_.structure[trackId] = landmarks[trackId];
                }
            }
            // Save outlier residual information
            Histogram<double> histoResiduals;
            std::cout << "\n"
                      << "=========================\n"
                      << " MSE Residual InitialPair Inlier:\n";
            ComputeResidualsHistogram(&histoResiduals);
            std::cout << "=========================" << std::endl;

            if (!sLogging_file_.empty())
            {
                using namespace htmlDocument;
                html_doc_stream_->pushInfo(htmlMarkup("h1","Essential Matrix."));
                std::ostringstream os;
                os << std::endl
                   << "-------------------------------" << "<br>"
                   << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
                   << view_I->s_Img_path << ","
                   << view_J->s_Img_path << "<br>"
                   << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
                   << "-- Resection status: " << "OK" << "<br>"
                   << "-- Nb points used for robust Essential matrix estimation: "
                   << xI.cols() << "<br>"
                   << "-- Nb points validated by robust estimation: "
                   << sfm_data_.structure.size() << "<br>"
                   << "-- % points validated: "
                   << sfm_data_.structure.size()/static_cast<float>(xI.cols())
                   << "<br>"
                   << "-------------------------------" << "<br>";
                html_doc_stream_->pushInfo(os.str());

                html_doc_stream_->pushInfo(htmlMarkup("h2",
                                                      "Residual of the robust estimation (Initial triangulation). Thresholded at: "
                                                      + toString(relativePose_info.found_residual_precision)));

                html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of residuals"));

                const std::vector<double> xBin = histoResiduals.GetXbinsValue();
                const auto range = autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

                htmlDocument::JSXGraphWrapper jsxGraph;
                jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
                jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
                jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
                                 relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
                jsxGraph.UnsuspendUpdate();
                jsxGraph.setViewport(range);
                jsxGraph.close();
                html_doc_stream_->pushInfo(jsxGraph.toStr());

                html_doc_stream_->pushInfo("<hr>");

                std::ofstream htmlFileStream( std::string(stlplus::folder_append_separator(sOut_directory_) +
                                                          "Reconstruction_Report.html").c_str());
                htmlFileStream << html_doc_stream_->getDoc();
            }
        }
        return !sfm_data_.structure.empty();
    }

    bool SequentialVISfMReconstructionEngine::AutomaticInitialPairChoice(Pair & initial_pair) const
    {
        // select a pair that have the largest baseline (mean angle between its bearing vectors).

        const unsigned iMin_inliers_count = 100;
        const float fRequired_min_angle = 3.0f;
        const float fLimit_max_angle = 60.0f; // More than 60 degree, we cannot rely on matches for initial pair seeding

        // List Views that support valid intrinsic (view that could be used for Essential matrix computation)
        std::set<IndexT> valid_views;
        for (Views::const_iterator it = sfm_data_.GetViews().begin();
             it != sfm_data_.GetViews().end(); ++it)
        {
            const View * v = it->second.get();
            if (sfm_data_.GetIntrinsics().count(v->id_intrinsic))
                valid_views.insert(v->id_view);
        }

        if (valid_views.size() < 2)
        {
            return false; // There is not view that support valid intrinsic data
        }

        std::vector<std::pair<double, Pair>> scoring_per_pair;

        // Compute the relative pose & the 'baseline score'
        C_Progress_display my_progress_bar( matches_provider_->pairWise_matches_.size(),
                                            std::cout,
                                            "Automatic selection of an initial pair:\n" );
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
        for (const std::pair<Pair, IndMatches> & match_pair : matches_provider_->pairWise_matches_)
        {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
            {
                ++my_progress_bar;

                const Pair current_pair = match_pair.first;

                const uint32_t I = std::min(current_pair.first, current_pair.second);
                const uint32_t J = std::max(current_pair.first, current_pair.second);
                if( (J - I) < 40 )  // TODO xinli
                if (valid_views.count(I) && valid_views.count(J))
                {
                    const View
                            * view_I = sfm_data_.GetViews().at(I).get(),
                            * view_J = sfm_data_.GetViews().at(J).get();
                    const Intrinsics::const_iterator
                            iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
                            iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

                    const auto
                            cam_I = iterIntrinsic_I->second.get(),
                            cam_J = iterIntrinsic_J->second.get();
                    if (cam_I && cam_J)
                    {
                        openMVG::tracks::STLMAPTracks map_tracksCommon;
                        shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

                        // Copy points correspondences to arrays for relative pose estimation
                        const size_t n = map_tracksCommon.size();
                        Mat xI(2,n), xJ(2,n);
                        size_t cptIndex = 0;
                        for (const auto & track_iter : map_tracksCommon)
                        {
                            auto iter = track_iter.second.cbegin();
                            const uint32_t i = iter->second;
                            const uint32_t j = (++iter)->second;

                            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
                            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
                            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
                            xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
                            ++cptIndex;
                        }

                        // Robust estimation of the relative pose
                        RelativePose_Info relativePose_info;
                        relativePose_info.initial_residual_tolerance = Square(4.0);

                        if (robustRelativePose(
                                cam_I, cam_J,
                                xI, xJ, relativePose_info,
                                {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
                                256)
                            && relativePose_info.vec_inliers.size() > iMin_inliers_count)
                        {
                            // Triangulate inliers & compute angle between bearing vectors
                            std::vector<float> vec_angles;
                            vec_angles.reserve(relativePose_info.vec_inliers.size());
                            const Pose3 pose_I = Pose3(Mat3::Identity(), Vec3::Zero());
                            const Pose3 pose_J = relativePose_info.relativePose;
                            for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
                            {
                                openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
                                std::advance(iterT, inlier_idx);
                                tracks::submapTrack::const_iterator iter = iterT->second.begin();
                                const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
                                const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
                                vec_angles.push_back(AngleBetweenRay(pose_I, cam_I, pose_J, cam_J,
                                                                     cam_I->get_ud_pixel(featI), cam_J->get_ud_pixel(featJ)));
                            }
                            // Compute the median triangulation angle
                            const unsigned median_index = vec_angles.size() / 2;
                            std::nth_element(
                                    vec_angles.begin(),
                                    vec_angles.begin() + median_index,
                                    vec_angles.end());
                            const float scoring_angle = vec_angles[median_index];
                            // Store the pair iff the pair is in the asked angle range [fRequired_min_angle;fLimit_max_angle]
                            if (scoring_angle > fRequired_min_angle &&
                                scoring_angle < fLimit_max_angle)
                            {
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
                                scoring_per_pair.emplace_back(scoring_angle, current_pair);
                            }
                        }
                    }
                }
            } // omp section
        }
        std::sort(scoring_per_pair.begin(), scoring_per_pair.end());
        // Since scoring is ordered in increasing order, reverse the order
        std::reverse(scoring_per_pair.begin(), scoring_per_pair.end());
        if (!scoring_per_pair.empty())
        {
            initial_pair = scoring_per_pair.begin()->second;
            return true;
        }
        return false;
    }

    /// Select a candidate initial pair
    bool SequentialVISfMReconstructionEngine::ChooseInitialPair(Pair & initialPairIndex) const
    {
        if (initial_pair_ != Pair(0,0))
        {
            // Internal initial pair is already initialized (so return it)
            initialPairIndex = initial_pair_;
        }
        else
        {
            // List Views that supports valid intrinsic
            std::set<IndexT> valid_views;
            for (const auto & view : sfm_data_.GetViews())
            {
                const View * v = view.second.get();
                if (sfm_data_.GetIntrinsics().find(v->id_intrinsic) != sfm_data_.GetIntrinsics().end())
                    valid_views.insert(v->id_view);
            }

            if (sfm_data_.GetIntrinsics().empty() || valid_views.empty())
            {
                std::cerr
                        << "There is no defined intrinsic data in order to compute an essential matrix for the initial pair."
                        << std::endl;
                return false;
            }

            std::cout << std::endl
                      << "----------------------------------------------------\n"
                      << "SequentialSfMReconstructionEngine::ChooseInitialPair\n"
                      << "----------------------------------------------------\n"
                      << " Pairs that have valid intrinsic and high support of points are displayed:\n"
                      << " Choose one pair manually by typing the two integer indexes\n"
                      << "----------------------------------------------------\n"
                      << std::endl;

            // Try to list the 10 top pairs that have:
            //  - valid intrinsics,
            //  - valid estimated Fundamental matrix.
            std::vector<uint32_t > vec_NbMatchesPerPair;
            std::vector<openMVG::matching::PairWiseMatches::const_iterator> vec_MatchesIterator;
            const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
            for (openMVG::matching::PairWiseMatches::const_iterator
                         iter = map_Matches.begin();
                 iter != map_Matches.end(); ++iter)
            {
                const Pair current_pair = iter->first;
                if (valid_views.count(current_pair.first) &&
                    valid_views.count(current_pair.second) )
                {
                    vec_NbMatchesPerPair.push_back(iter->second.size());
                    vec_MatchesIterator.push_back(iter);
                }
            }
            // sort the Pairs in descending order according their correspondences count
            using namespace stl::indexed_sort;
            std::vector<sort_index_packet_descend<uint32_t, uint32_t>> packet_vec(vec_NbMatchesPerPair.size());
            sort_index_helper(packet_vec, &vec_NbMatchesPerPair[0], std::min((size_t)10, vec_NbMatchesPerPair.size()));

            for (size_t i = 0; i < std::min((size_t)10, vec_NbMatchesPerPair.size()); ++i) {
                const uint32_t index = packet_vec[i].index;
                openMVG::matching::PairWiseMatches::const_iterator iter = vec_MatchesIterator[index];
                std::cout << "(" << iter->first.first << "," << iter->first.second <<")\t\t"
                          << iter->second.size() << " matches" << std::endl;
            }

            // Ask the user to choose an initial pair (by set some view ids)
            std::cout << std::endl << " type INITIAL pair ids: X enter Y enter\n";
            int val, val2;
            if ( std::cin >> val && std::cin >> val2) {
                initialPairIndex.first = val;
                initialPairIndex.second = val2;
            }
        }

        std::cout << "\nPutative starting pair is: (" << initialPairIndex.first
                  << "," << initialPairIndex.second << ")" << std::endl;

        // Check validity of the initial pair indices:
        if (features_provider_->feats_per_view.find(initialPairIndex.first) == features_provider_->feats_per_view.end() ||
            features_provider_->feats_per_view.find(initialPairIndex.second) == features_provider_->feats_per_view.end())
        {
            std::cerr << "At least one of the initial pair indices is invalid."
                      << std::endl;
            return false;
        }
        return true;
    }

    /// Functor to sort a vector of pair given the pair's second value
    template<class T1, class T2, class Pred = std::less<T2>>
    struct sort_pair_second {
        bool operator()(const std::pair<T1,T2>&left,
                        const std::pair<T1,T2>&right)
        {
            Pred p;
            return p(left.second, right.second);
        }
    };

    /**
     * @brief Estimate images on which we can compute the resectioning safely.
     *
     * @param[out] vec_possible_indexes: list of indexes we can use for resectioning.
     * @return False if there is no possible resection.
     *
     * Sort the images by the number of features id shared with the reconstruction.
     * Select the image I that share the most of correspondences.
     * Then keep all the images that have at least:
     *  0.75 * #correspondences(I) common correspondences to the reconstruction.
     */
    bool SequentialVISfMReconstructionEngine::FindImagesWithPossibleResection(
            std::vector<uint32_t> & vec_possible_indexes)
    {
        // Threshold used to select the best images
        static const float dThresholdGroup = 0.75f;

        // TODO xinDEBUG
//        std::cout <<"set_remaining_view_id_ contains : \n";
//        for( auto&view_id:set_remaining_view_id_ )
//        {
//            std::cout << view_id << " ";
//        }
//        std::cout << std::endl;
//
//        std::cout <<"sfm_data_.views first contains : \n";
//        for( auto&view_id:sfm_data_.views )
//        {
//            std::cout << view_id.first << " ";
//        }
//        std::cout << std::endl;
//        std::cout <<"sfm_data_.views second viewId contains : \n";
//        for( auto&view_id:sfm_data_.views )
//        {
//            std::cout << view_id.second->id_view << " ";
//        }
//        std::cout << std::endl;

        vec_possible_indexes.clear();

        // check index view == index pose
        IndexT left_view_id, right_view_id;
        IndexT left_pose_id, right_pose_id;
        {
            auto it = sfm_data_.views.begin();
            auto c_it = sfm_data_.views.rbegin();
            left_view_id = it->first;
            right_view_id = c_it->first;
        }

        std::cout << "left_view_id = " << left_view_id << std::endl;
        std::cout << "right_view_id = " << right_view_id << std::endl;

        // TODO xinli make sure pose_id == view_id
        {
            auto it = sfm_data_.poses.begin();
            auto c_it = sfm_data_.poses.rbegin();
            left_pose_id = it->first;
            right_pose_id = c_it->first;

            std::cout << "left_pose_id = " << left_pose_id << std::endl;
            std::cout << "right_pose_id = " << right_pose_id << std::endl;
            std::cout << "sfm_data_.poses.size() = " << sfm_data_.poses.size() << std::endl;

            if( (left_pose_id - left_view_id) < 10 ) left_pose_id = left_view_id;
            else left_pose_id -= 10;
            if( (right_view_id - right_pose_id) < 10 ) right_pose_id = right_view_id;
            else right_pose_id += 10;
        }

        std::cout << "left_pose_id = " << left_pose_id << std::endl;
        std::cout << "right_pose_id = " << right_pose_id << std::endl;

        if (set_remaining_view_id_.empty() || sfm_data_.GetLandmarks().empty())
            return false;

        // Collect tracksIds
        std::set<uint32_t> reconstructed_trackId;
        std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
                       std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
                       stl::RetrieveKey());

        Pair_Vec vec_putative; // ImageId, NbPutativeCommonPoint
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
        for (std::set<uint32_t>::const_iterator iter = set_remaining_view_id_.begin();
             iter != set_remaining_view_id_.end(); ++iter)
        {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
            {
                const uint32_t viewId = *iter;

                if( viewId >= left_pose_id && viewId <= right_pose_id )
                {
                    // Compute 2D - 3D possible content
                    openMVG::tracks::STLMAPTracks map_tracksCommon;
                    shared_track_visibility_helper_->GetTracksInImages({viewId}, map_tracksCommon);

                    // TODO xinDEBUG
//                    std::cout << "viewId = " << viewId << " map_tracksCommon.size() : " << map_tracksCommon.size() << std::endl;

                    if (!map_tracksCommon.empty())
                    {
                        std::set<uint32_t> set_tracksIds;
                        tracks::TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

                        // Count the common possible putative point
                        //  with the already 3D reconstructed trackId
                        std::vector<uint32_t> vec_trackIdForResection;
                        std::set_intersection(set_tracksIds.cbegin(), set_tracksIds.cend(),
                                              reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
                                              std::back_inserter(vec_trackIdForResection));

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
                        {
                            vec_putative.emplace_back(viewId, vec_trackIdForResection.size());
                        }
                    }
                }


            }
        }

        // Sort by the number of matches to the 3D scene.
        std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<uint32_t, uint32_t, std::greater<uint32_t>>());


        // TODO xinDEBUG
//        std::cout <<"vec_putative contains : ";
//        for( auto&view_id:vec_putative )
//        {
//            std::cout << view_id.second<< " ";
//        }
//        std::cout << std::endl;

        // If the list is empty or if the list contains images with no correspdences
        // -> (no resection will be possible)
        if (vec_putative.empty() || vec_putative[0].second == 0)
        {
            // All remaining images cannot be used for pose estimation
            // TODO xinDEBUG
            std::cout <<"All remaining images cannot be used for pose estimation" << std::endl;
            set_remaining_view_id_.clear();
            return false;
        }

        // Add the image view index that share the most of 2D-3D correspondences
        vec_possible_indexes.push_back(vec_putative[0].first);

        // Then, add all the image view indexes that have at least N% of the number of the matches of the best image.
        const IndexT M = vec_putative[0].second; // Number of 2D-3D correspondences
        const size_t threshold = static_cast<uint32_t>(dThresholdGroup * M);
        for (size_t i = 1; i < vec_putative.size() &&
                           vec_putative[i].second > threshold; ++i)
        {
            vec_possible_indexes.push_back(vec_putative[i].first);
        }
        return true;
    }
     /**
     * @brief Estimate images on which we can compute the resectioning safely.
     *
     * @param[out] vec_possible_indexes: list of indexes we can use for resectioning.
     * @return False if there is no possible resection.
     *
     * Sort the images by the number of features id shared with the reconstruction.
     * Select the image I that share the most of correspondences.
     * Then keep all the images that have at least:
     *  0.75 * #correspondences(I) common correspondences to the reconstruction.
     */
    bool SequentialVISfMReconstructionEngine::FindImagesWithPossibleResection_VIinit(
            std::vector<uint32_t> & vec_possible_indexes)
    {
        // Threshold used to select the best images
        static const float dThresholdGroup = 0.75f;

        vec_possible_indexes.clear();

        if (set_remaining_view_id_vi_init_.empty() || sfm_data_.GetLandmarks().empty())
            return false;

        // Collect tracksIds
        std::set<uint32_t> reconstructed_trackId;
        std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
                       std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
                       stl::RetrieveKey());

        Pair_Vec vec_putative; // ImageId, NbPutativeCommonPoint
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
        for (std::set<uint32_t>::const_iterator iter = set_remaining_view_id_vi_init_.begin();
             iter != set_remaining_view_id_vi_init_.end(); ++iter)
        {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
            {
                const uint32_t viewId = *iter;

                // Compute 2D - 3D possible content
                openMVG::tracks::STLMAPTracks map_tracksCommon;
                shared_track_visibility_helper_->GetTracksInImages({viewId}, map_tracksCommon);

                if (!map_tracksCommon.empty())
                {
                    std::set<uint32_t> set_tracksIds;
                    tracks::TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

                    // Count the common possible putative point
                    //  with the already 3D reconstructed trackId
                    std::vector<uint32_t> vec_trackIdForResection;
                    std::set_intersection(set_tracksIds.cbegin(), set_tracksIds.cend(),
                                          reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
                                          std::back_inserter(vec_trackIdForResection));

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
                    {
                        vec_putative.emplace_back(viewId, vec_trackIdForResection.size());
                    }
                }
            }
        }

        // Sort by the number of matches to the 3D scene.
        std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<uint32_t, uint32_t, std::greater<uint32_t>>());

        // If the list is empty or if the list contains images with no correspdences
        // -> (no resection will be possible)
        if (vec_putative.empty() || vec_putative[0].second == 0)
        {
            // All remaining images cannot be used for pose estimation
            set_remaining_view_id_vi_init_.clear();
            return false;
        }

        // Add the image view index that share the most of 2D-3D correspondences
        vec_possible_indexes.push_back(vec_putative[0].first);

        // Then, add all the image view indexes that have at least N% of the number of the matches of the best image.
        const IndexT M = vec_putative[0].second; // Number of 2D-3D correspondences
        const size_t threshold = static_cast<uint32_t>(dThresholdGroup * M);
        for (size_t i = 1; i < vec_putative.size() &&
                           vec_putative[i].second > threshold; ++i)
        {
            vec_possible_indexes.push_back(vec_putative[i].first);
        }
        return true;
    }

    /**
     * @brief Add one image to the 3D reconstruction. To the resectioning of
     * the camera and triangulate all the new possible tracks.
     * @param[in] viewIndex: image index to add to the reconstruction.
     *
     * A. Compute 2D/3D matches
     * B. Look if intrinsic data is known or not
     * C. Do the resectioning: compute the camera pose.
     * D. Refine the pose of the found camera
     * E. Update the global scene with the new camera
     * F. Update the observations into the global scene structure
     * G. Triangulate new possible 2D tracks
     */
    bool SequentialVISfMReconstructionEngine::Resection(const uint32_t viewIndex)
    {
        using namespace tracks;

        // A. Compute 2D/3D matches
        // A1. list tracks ids used by the view
        openMVG::tracks::STLMAPTracks map_tracksCommon;
        shared_track_visibility_helper_->GetTracksInImages({viewIndex}, map_tracksCommon);
        std::set<uint32_t> set_tracksIds;
        TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

        // A2. intersects the track list with the reconstructed
        std::set<uint32_t> reconstructed_trackId;
        std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
                       std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
                       stl::RetrieveKey());

        // Get the ids of the already reconstructed tracks
        std::set<uint32_t> set_trackIdForResection;
        std::set_intersection(set_tracksIds.cbegin(), set_tracksIds.cend(),
                              reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
                              std::inserter(set_trackIdForResection, set_trackIdForResection.begin()));

        if (set_trackIdForResection.empty())
        {
            // No match. The image has no connection with already reconstructed points.
            std::cout << std::endl
                      << "-------------------------------" << "\n"
                      << "-- Resection of camera index: " << viewIndex << "\n"
                      << "-- Resection status: " << "FAILED" << "\n"
                      << "-------------------------------" << std::endl;
            return false;
        }

        // Get back featId associated to a tracksID already reconstructed.
        // These 2D/3D associations will be used for the resection.
        std::vector<uint32_t> vec_featIdForResection;
        TracksUtilsMap::GetFeatIndexPerViewAndTrackId(map_tracksCommon,
                                                      set_trackIdForResection,
                                                      viewIndex,
                                                      &vec_featIdForResection);

        // Localize the image inside the SfM reconstruction
        Image_Localizer_Match_Data resection_data;
        resection_data.pt2D.resize(2, set_trackIdForResection.size());
        resection_data.pt3D.resize(3, set_trackIdForResection.size());

        // B. Look if the intrinsic data is known or not
        const View * view_I = sfm_data_.GetViews().at(viewIndex).get();
        std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic;
        if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic))
        {
            optional_intrinsic = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic);
        }

        // Setup the track 2d observation for this new view
        Mat2X pt2D_original(2, set_trackIdForResection.size());
        std::set<uint32_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
        std::vector<uint32_t>::const_iterator iterfeatId = vec_featIdForResection.begin();
        for (size_t cpt = 0; cpt < vec_featIdForResection.size(); ++cpt, ++iterTrackId, ++iterfeatId)
        {
            resection_data.pt3D.col(cpt) = sfm_data_.GetLandmarks().at(*iterTrackId).X;
            resection_data.pt2D.col(cpt) = pt2D_original.col(cpt) =
                    features_provider_->feats_per_view.at(viewIndex)[*iterfeatId].coords().cast<double>();
            // Handle image distortion if intrinsic is known (to ease the resection)
            if (optional_intrinsic && optional_intrinsic->have_disto())
            {
                resection_data.pt2D.col(cpt) = optional_intrinsic->get_ud_pixel(resection_data.pt2D.col(cpt));
            }
        }

        // C. Do the resectioning: compute the camera pose
        std::cout << std::endl
                  << "-------------------------------" << std::endl
                  << "-- Robust Resection of view: " << viewIndex << std::endl;

        geometry::Pose3 pose;
        const bool bResection = sfm::SfM_Localizer::Localize
                (
                        optional_intrinsic ? resection::SolverType::P3P_NORDBERG_ECCV18 : resection::SolverType::DLT_6POINTS,
                        {view_I->ui_width, view_I->ui_height},
                        optional_intrinsic.get(),
                        resection_data,
                        pose
                );
        resection_data.pt2D = std::move(pt2D_original); // restore original image domain points

        if (!sLogging_file_.empty())
        {
            using namespace htmlDocument;
            std::ostringstream os;
            os << "Resection of Image index: <" << viewIndex << "> image: "
               << view_I->s_Img_path <<"<br> \n";
            html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

            os.str("");
            os << std::endl
               << "-------------------------------" << "<br>"
               << "-- Robust Resection of camera index: <" << viewIndex << "> image: "
               <<  view_I->s_Img_path <<"<br>"
               << "-- Threshold: " << resection_data.error_max << "<br>"
               << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "<br>"
               << "-- Nb points used for Resection: " << vec_featIdForResection.size() << "<br>"
               << "-- Nb points validated by robust estimation: " << resection_data.vec_inliers.size() << "<br>"
               << "-- % points validated: "
               << resection_data.vec_inliers.size()/static_cast<float>(vec_featIdForResection.size()) << "<br>"
               << "-------------------------------" << "<br>";
            html_doc_stream_->pushInfo(os.str());
        }

        if (!bResection)
            return false;

        // D. Refine the pose of the found camera.
        // We use a local scene with only the 3D points and the new camera.
        {
            const bool b_new_intrinsic = (optional_intrinsic == nullptr);
            // A valid pose has been found (try to refine it):
            // If no valid intrinsic as input:
            //  init a new one from the projection matrix decomposition
            // Else use the existing one and consider it as constant.
            if (b_new_intrinsic)
            {
                // setup a default camera model from the found projection matrix
                Mat3 K, R;
                Vec3 t;
                KRt_From_P(resection_data.projection_matrix, &K, &R, &t);

                const double focal = (K(0,0) + K(1,1))/2.0;
                const Vec2 principal_point(K(0,2), K(1,2));

                // Create the new camera intrinsic group
                switch (cam_type_)
                {
                    case PINHOLE_CAMERA:
                        optional_intrinsic =
                                std::make_shared<Pinhole_Intrinsic>
                                        (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
                        break;
                    case PINHOLE_CAMERA_RADIAL1:
                        optional_intrinsic =
                                std::make_shared<Pinhole_Intrinsic_Radial_K1>
                                        (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
                        break;
                    case PINHOLE_CAMERA_RADIAL3:
                        optional_intrinsic =
                                std::make_shared<Pinhole_Intrinsic_Radial_K3>
                                        (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
                        break;
                    case PINHOLE_CAMERA_BROWN:
                        optional_intrinsic =
                                std::make_shared<Pinhole_Intrinsic_Brown_T2>
                                        (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
                        break;
                    case PINHOLE_CAMERA_FISHEYE:
                        optional_intrinsic =
                                std::make_shared<Pinhole_Intrinsic_Fisheye>
                                        (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
                        break;
                    default:
                        std::cerr << "Try to create an unknown camera type." << std::endl;
                        return false;
                }
            }
            const bool b_refine_pose = true;
            const bool b_refine_intrinsics = false;
            if (!sfm::SfM_Localizer::RefinePose(
                    optional_intrinsic.get(), pose,
                    resection_data, b_refine_pose, b_refine_intrinsics))
            {
                return false;
            }

            // E. Update the global scene with:
            // - the new found camera pose
            sfm_data_.poses[view_I->id_pose] = pose;
            // - track the view's AContrario robust estimation found threshold
            map_ACThreshold_.insert({viewIndex, resection_data.error_max});
            // - intrinsic parameters (if the view has no intrinsic group add a new one)
            if (b_new_intrinsic)
            {
                // Since the view have not yet an intrinsic group before, create a new one
                IndexT new_intrinsic_id = 0;
                if (!sfm_data_.GetIntrinsics().empty())
                {
                    // Since some intrinsic Id already exists,
                    //  we have to create a new unique identifier following the existing one
                    std::set<IndexT> existing_intrinsicId;
                    std::transform(sfm_data_.GetIntrinsics().cbegin(), sfm_data_.GetIntrinsics().cend(),
                                   std::inserter(existing_intrinsicId, existing_intrinsicId.begin()),
                                   stl::RetrieveKey());
                    new_intrinsic_id = (*existing_intrinsicId.rbegin())+1;
                }
                sfm_data_.views.at(viewIndex)->id_intrinsic = new_intrinsic_id;
                sfm_data_.intrinsics[new_intrinsic_id] = optional_intrinsic;
            }
        }

        // F. List tracks that share content with this view and add observations and new 3D track if required.
        //    - If the track already exists (look if the new view tracks observation are valid)
        //    - If the track does not exists, try robust triangulation & add the new valid view track observation
        {
            // Get information of new view
            const IndexT I = viewIndex;
            const View * view_I = sfm_data_.GetViews().at(I).get();
            const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
            const Pose3 pose_I = sfm_data_.GetPoseOrDie(view_I);

            // Vector of all already reconstructed views
            const std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);

            // Go through each track and look if we must add new view observations or new 3D points
            for (const std::pair<uint32_t, tracks::submapTrack>& trackIt : map_tracksCommon)
            {
                const uint32_t trackId = trackIt.first;
                const tracks::submapTrack & track = trackIt.second;

                // List the potential view observations of the track
                const tracks::submapTrack & allViews_of_track = map_tracks_[trackId];

                // List to save the new view observations that must be added to the track
                std::set<IndexT> new_track_observations_valid_views;

                // If the track was already reconstructed
                if (sfm_data_.structure.count(trackId) != 0)
                {
                    // Since the 3D point was triangulated before we add the new the Inth view observation
                    new_track_observations_valid_views.insert(I);
                }
                else
                {
                    // Go through the views that observe this track & look if a successful triangulation can be done
                    for (const std::pair<IndexT, IndexT>& trackViewIt : allViews_of_track)
                    {
                        const IndexT & J = trackViewIt.first;
                        // If view is valid try triangulation
                        if (J != I && valid_views.count(J) != 0)
                        {
                            // If successfuly triangulated add the observation from J view
                            if (sfm_data_.structure.count(trackId) != 0)
                            {
                                new_track_observations_valid_views.insert(J);
                            }
                            else
                            {
                                const View * view_J = sfm_data_.GetViews().at(J).get();
                                const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
                                const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
                                const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();

                                // Position of the point in view I
                                const Vec2 xI = features_provider_->feats_per_view.at(I)[track.at(I)].coords().cast<double>();

                                // Try to triangulate a 3D point from J view
                                // A new 3D point must be added
                                // Triangulate it
                                const Vec2
                                        xI_ud = cam_I->get_ud_pixel(xI),
                                        xJ_ud = cam_J->get_ud_pixel(xJ);
                                Vec3 X = Vec3::Zero();

                                if (Triangulate2View(
                                        pose_I.rotation(),
                                        pose_I.translation(),
                                        (*cam_I)(xI_ud),
                                        pose_J.rotation(),
                                        pose_J.translation(),
                                        (*cam_J)(xJ_ud),
                                        X,
                                        triangulation_method_))
                                {
                                    // Check triangulation result
                                    const double angle = AngleBetweenRay(
                                            pose_I, cam_I, pose_J, cam_J, xI_ud, xJ_ud);
                                    const Vec2 residual_I = cam_I->residual(pose_I(X), xI);
                                    const Vec2 residual_J = cam_J->residual(pose_J(X), xJ);
                                    if (
                                        //  - Check angle (small angle leads to imprecise triangulation)
                                            angle > 2.0 &&
                                            //  - Check residual values (must be inferior to the found view's AContrario threshold)
                                            residual_I.norm() < std::max(4.0, map_ACThreshold_.at(I)) &&
                                            residual_J.norm() < std::max(4.0, map_ACThreshold_.at(J))
                                        // Cheirality as been tested already in Triangulate2View
                                            )
                                    {
                                        // Add a new track
                                        Landmark & landmark = sfm_data_.structure[trackId];
                                        landmark.X = X;
                                        new_track_observations_valid_views.insert(I);
                                        new_track_observations_valid_views.insert(J);
                                    } // 3D point is valid
                                }
                                else
                                {
                                    // We mark the view to add the observations once the point is triangulated
                                    new_track_observations_valid_views.insert(J);
                                } // 3D point is invalid
                            }
                        }
                    }// Go through all the views
                }// If new point

                // If successfuly triangulated, add the valid view observations
                if (sfm_data_.structure.count(trackId) != 0 &&
                    !new_track_observations_valid_views.empty())
                {
                    Landmark & landmark = sfm_data_.structure[trackId];
                    // Check if view feature point observations of the track are valid (residual, depth) or not
                    for (const IndexT & J: new_track_observations_valid_views)
                    {
                        const View * view_J = sfm_data_.GetViews().at(J).get();
                        const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
                        const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
                        const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();
                        const Vec2 xJ_ud = cam_J->get_ud_pixel(xJ);

                        const Vec2 residual = cam_J->residual(pose_J(landmark.X), xJ);
                        if (CheiralityTest((*cam_J)(xJ_ud), pose_J, landmark.X)
                            && residual.norm() < std::max(4.0, map_ACThreshold_.at(J))
                                )
                        {
                            landmark.obs[J] = Observation(xJ, allViews_of_track.at(J));
                        }
                    }
                }
            }// All the tracks in the view
        }
        return true;
    }
}
}