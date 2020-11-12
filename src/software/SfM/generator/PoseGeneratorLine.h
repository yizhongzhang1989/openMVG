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

    PoseGeneratorLine(int deltaT, int deltaT_IMU, bool storeIMU = false)
    : t_cam_ms(0), t_imu_ms(0), storeIMU_(storeIMU), deltaT(deltaT), deltaT_IMU(deltaT_IMU)
    {
        IMUs.clear();
    }
    Pose Generate() override
    {
        double t_cam = 1e-3 * t_cam_ms;
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
            while(t_imu_ms <= t_cam_ms)
            {
                double t_imu = 1e-3 * t_imu_ms;
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
            }
        }

        t_cam_ms += deltaT;

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

    int getDeltaT() const override
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
    // sampling period in ms
    int deltaT;
    int deltaT_IMU;
    // current time in ms
    int t_cam_ms;
    int t_imu_ms;
    // IMU measurements
    IMUMeasurements IMUs;
    bool storeIMU_;
};

}  // namespace generator


#endif  // POSE_GENERATOR_LINE_H_