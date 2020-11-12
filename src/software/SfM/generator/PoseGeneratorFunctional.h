#pragma once
#ifndef POSE_GENERATOR_FUNCTIONAL_H_
#define POSE_GENERATOR_FUNCTIONAL_H_

#include <functional>
#include "PoseGeneratorBase.h"
#include "types.h"

namespace generator
{

class PoseGeneratorFunctional : public PoseGeneratorBase<Pose, Eigen::aligned_allocator<Pose>>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef std::function<InversePose (double t)> PoseFunction;

    PoseGeneratorFunctional(int deltaT, int deltaT_IMU, PoseFunction& func, bool storeIMU = false)
            : t_cam_ms(0), t_imu_ms(0), storeIMU_(storeIMU), deltaT(deltaT), deltaT_IMU(deltaT_IMU), poseFunction(func)
    {
        IMUs.clear();
    }

    Pose Generate() override
    {
        static const double eps = 1e-6;

        double t_cam = 1e-3 * t_cam_ms;
        InversePose inv_pose = poseFunction(t_cam);
        Pose pose;
        pose.q = inv_pose.q.inverse();
        pose.t = - (inv_pose.q.inverse() * inv_pose.p);

        while(t_imu_ms <= t_cam_ms)
        {
            double t_imu = 1e-3 * t_imu_ms;

            // differentiate pose
            InversePose p_t = poseFunction(t_imu);
            InversePose p_t_1 = poseFunction(t_imu + eps);
            InversePose p_t_2 = poseFunction(t_imu + 2 * eps);

            Eigen::Vector3d v_t = (p_t_1.p - p_t.p) / eps;
            Eigen::Vector3d v_t_1 = (p_t_2.p - p_t_1.p) / eps;
            Eigen::Vector3d a_t = (v_t_1 - v_t) / eps;

            Eigen::Quaterniond dq = p_t.q.inverse() * p_t_1.q;
            Eigen::Vector3d theta = 2 * dq.vec();
            Eigen::Vector3d omega = - theta / eps;

            Eigen::Vector3d a_local = pose.q * a_t;
            IMUs.emplace_back(a_local,omega,t_imu_ms);

            t_imu_ms + deltaT_IMU;
        }

        t_cam_ms + deltaT;

        return pose;
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
    PoseFunction& poseFunction;
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

#endif  // POSE_GENERATOR_FUNCTIONAL_H_