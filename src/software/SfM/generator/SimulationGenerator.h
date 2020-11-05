#pragma once
#ifndef SIMULATION_GENERATOR_H_
#define SIMULATION_GENERATOR_H_

#include <utility>

#include "PointGenerator.h"
#include "PoseGeneratorBase.h"
#include "CameraBase.h"
#include "ColorGenerator.h"
#include "types.h"

namespace generator
{

class SimulationGenerator
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    struct SimulationConfig
    {
        int n_points;
        int n_lines;
        int n_poses;
        int image_width;
        int image_height;
        SimulationConfig()
        {
            n_points = 2000;
            n_lines = 0;
            n_poses = 100;
            image_width = 1024;
            image_height = 768;
        }
    };
    struct NoiseConfig
    {
        double stddev_rotation;
        double stddev_translation;
        double stddev_point3d;
        double stddev_point2d;
        bool perturb_rotation;
        bool perturb_translation;
        bool perturb_point3d;
        bool perturb_point2d;
        NoiseConfig()
        {
            stddev_rotation = 0.001;
            stddev_translation = 0.01;
            stddev_point3d = 0.01;
            stddev_point2d = 1.0;
            perturb_rotation = true;
            perturb_translation = true;
            perturb_point3d = false;
            perturb_point2d = true;
        }
    };
    enum PointGenerationMode
    {
        RANDOM_POINT,
        SAMPLING_SURFACE
    };
    SimulationGenerator(PoseGeneratorBase<Pose>* pPoseGenerator,
                        PointGenerator* pPointGenerator,
                        CameraBase<Eigen::Vector3d,Eigen::Vector2d>* pCamera)
    :mpPoseGenerator(pPoseGenerator), mpPointGenerator(pPointGenerator),
    mpCamera(pCamera), mPointMode(RANDOM_POINT)
    {
        ;
    }
    SimulationGenerator(PoseGeneratorBase<Pose>* pPoseGenerator,
                        std::string ObjFileName,
                        CameraBase<Eigen::Vector3d,Eigen::Vector2d>* pCamera)
            :mpPoseGenerator(pPoseGenerator), strObjFileName(std::move(ObjFileName)),
            mpCamera(pCamera), mPointMode(SAMPLING_SURFACE), mpPointGenerator(nullptr)
    {
        ;
    }
    bool Generate(Simulation_Data& sfm_data, const SimulationConfig& cfg);
    static void AddNoise(const Simulation_Data& sfm_data, Simulation_Data& sfm_data_noisy, NoiseConfig& cfg_noise);
    void Save(Simulation_Data& sfm_data, const std::string& outPath);
    void SaveIMU(const IMUMeasurements& imu_data, const std::string& imu_file);

    template<typename PoseGeneratorType>
    const IMUMeasurements& getIMUMeasurements() const
    {
        PoseGeneratorType* pPoseGenerator = dynamic_cast<PoseGeneratorType*>(mpPoseGenerator);
        if(pPoseGenerator && pPoseGenerator->hasIMU())
        {
            return pPoseGenerator->getIMUMeasurements();
        }
        else
            throw std::runtime_error("Can not get IMU measurements");
    }
private:
    PoseGeneratorBase<Pose>* mpPoseGenerator;
    PointGenerator* mpPointGenerator;
    CameraBase<Eigen::Vector3d,Eigen::Vector2d>* mpCamera;
    ColorGenerator<Color> mColorGenerator;
    const PointGenerationMode mPointMode;
    std::string strObjFileName;
};

}  // namespace generator

#endif  // SIMULATION_GENERATOR_H_