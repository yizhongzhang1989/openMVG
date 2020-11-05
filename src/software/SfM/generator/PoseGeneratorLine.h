#pragma once
#ifndef POSE_GENERATOR_LINE_H_
#define POSE_GENERATOR_LINE_H_

#include "PoseGeneratorBase.h"
#include "types.h"

namespace generator
{

class PoseGeneratorLine : public PoseGeneratorBase<Pose>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PoseGeneratorLine(double freq_img, double freq_imu, bool storeIMU)
    : t_cam(0.0), t_imu(0.0), storeIMU_(storeIMU)
    {
        IMUs.clear();

        deltaT = 1.0 / freq_img;
        deltaT_IMU = 1.0 / freq_imu;
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
                IMUs.emplace_back(acc,gyro,t_imu);

                t_imu += deltaT_IMU;
            }
        }

        t_cam += deltaT;

        return p;
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
    double deltaT;
    double deltaT_IMU;
    double t_cam;
    double t_imu;
    IMUMeasurements IMUs;
    bool storeIMU_;
};

}  // namespace generator


#endif  // POSE_GENERATOR_LINE_H_