#pragma once
#ifndef POSE_GENERATOR_CONST_ACC_H_
#define POSE_GENERATOR_CONST_ACC_H_

#include "PoseGeneratorBase.h"
#include "types.h"

namespace generator
{
class PoseGeneratorConstAcc : public PoseGeneratorBase<Pose>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    PoseGeneratorConstAcc(double acc_x, double acc_y, double acc_z, double freq_img, double freq_imu, bool storeIMU = false, LookDirection direction = FORWARD)
    : acc_x(acc_x), acc_y(acc_y), acc_z(acc_z), direction(direction), t_cam(0.0), t_imu(0.0), storeIMU_(storeIMU)
    {
        IMUs.clear();

        deltaT = 1.0 / freq_img;
        deltaT_IMU = 1.0 / freq_imu;
    }
    Pose Generate() override
    {
        static const Eigen::Vector3d gravity(0.0,0.0,-GRAVITY);
        static const double eps = 1e-6;

        Pose p;
        Eigen::Vector3d acc(acc_x,acc_y,acc_z);

        // generate position, t: local ->global, position of camera center, t_wc
        Eigen::Vector3d t = 0.5 * acc * t_cam * t_cam;

        // generate orientation, R: local -> global, R_wc
        Eigen::Matrix3d R;
        Eigen::Vector3d rx,ry,rz;

        Eigen::Vector3d left_vec;
        if(fabs(acc_x) > eps)
        {
            left_vec = Eigen::Vector3d(-acc_y/acc_x, 1.0, 0.0);
        }
        else if(fabs(acc_y) > eps)
        {
            left_vec = Eigen::Vector3d(1.0, -acc_x/acc_y, 0.0);
        }
        else
        {
            left_vec = Eigen::Vector3d(-t.x(),-t.y(),0.0);
        }

        switch(direction)
        {
            case FORWARD:
            {
                rz = acc;
                ry = left_vec;
                rx = ry.cross(rz);
            } break;
            case LEFTWARD:
            {
                rz = left_vec;
                ry = -acc;
                rx = ry.cross(rz);
            } break;
            default:
            {
                throw std::runtime_error("unrecognized direction.");
            }
        }

        rx = rx / rx.norm();
        ry = ry / ry.norm();
        rz = rz / rz.norm();

        R.col(0) = rx;
        R.col(1) = ry;
        R.col(2) = rz;

        // p.q: global -> local, q_cw
        p.q = Eigen::Quaterniond(R.transpose());
        // p.t: global -> local, t_cw
        p.t = - R.transpose() * t;

        Eigen::Vector3d acc_local = p.q * (acc + gravity);
        while(t_imu <= t_cam)
        {
            IMUs.emplace_back(acc_local, Eigen::Vector3d::Zero(), t_imu);
            t_imu += deltaT_IMU;
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
    double acc_x, acc_y, acc_z;
    double deltaT, deltaT_IMU;
    double t_cam, t_imu;
    IMUMeasurements IMUs;
    LookDirection direction;
    bool storeIMU_;
};
}

#endif  // POSE_GENERATOR_CONST_ACC_H_