#pragma once
#ifndef POSE_GENERATOR_LINE_H_
#define POSE_GENERATOR_LINE_H_

#include "PoseGeneratorBase.h"
#include "types.h"

namespace generator
{

class PoseGeneratorLine : public PoseGeneratorBase<Pose, Eigen::aligned_allocator<Pose>>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PoseGeneratorLine(double deltaT_s, int deltaT_IMU_ms, bool storeIMU = false)
    : t_cam(0.0), t_imu_ms(0), storeIMU_(storeIMU), deltaT(deltaT_s), deltaT_IMU(deltaT_IMU_ms)
    {
        IMUs.clear();
    }
    Pose Generate() override
    {
        double omega = 2 * PI * t_cam;
        double sin_omega = sin(omega);

        double x = 0.0;
        double y = 0.1 * sin_omega;
        double z = t_cam + 0.05 * sin_omega;
        Eigen::Vector3d t(x,y,z);

        Pose p;
        p.t = - t;
        p.q.setIdentity();

//        double cos_omega = cos(omega);
//        double dx = 0.0;
//        double dy = 0.1 * cos_omega * 2 * PI;
//        double dz = 1.0 + 0.05 * cos(omega) * 2 * PI;

        if(storeIMU_)
        {
            double t_imu = 1e-3 * t_imu_ms;
            while(t_imu <= t_cam)
            {
                double PI_2 = 2 * PI;
                double PI_2_2 = PI_2 * PI_2;
                double omega_imu = PI_2 * t_imu;
                double sin_omega_imu = sin(omega_imu);

                double dx2 = 0.0;
                double dy2 = - 0.1 * sin_omega_imu * PI_2_2;
                double dz2 = - 0.05 * sin_omega_imu * PI_2_2;
                Eigen::Vector3d acc(dx2,dy2,dz2);
                Eigen::Vector3d gyro(0.0,0.0,0.0);
                IMUs.emplace_back(acc,gyro,t_imu_ms);

                t_imu_ms += deltaT_IMU;
                t_imu = 1e-3 * t_imu_ms;
            }
        }

        t_cam += deltaT;

        return p;
    }

    STLVector<Pose> Generate(int num_poses) override
    {
        STLVector<Pose> poses;
        for(int i = 0; i < num_poses; i++)
        {
            poses.push_back(Generate());
        }
        return poses;
    }

    double getDeltaT() const override
    {
        return deltaT;
    }

    const IMUMeasurements& getIMUMeasurements() const
    {
        return IMUs;
    }

    bool hasIMU() const
    {
        return storeIMU_;
    }
private:
    // camera sampling period in s
    double deltaT;
    // IMU sampling period in ms
    int deltaT_IMU;
    // camera current time in s
    double t_cam;
    // IMU current time in ms
    int t_imu_ms;
    // IMU measurements
    IMUMeasurements IMUs;
    bool storeIMU_;
};

}  // namespace generator


#endif  // POSE_GENERATOR_LINE_H_