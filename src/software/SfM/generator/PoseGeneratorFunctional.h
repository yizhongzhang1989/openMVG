#pragma once
#ifndef POSE_GENERATOR_FUNCTIONAL_H_
#define POSE_GENERATOR_FUNCTIONAL_H_

#include <functional>
#include "PoseGeneratorBase.h"
#include "types.h"

namespace generator
{

class PoseGeneratorFunctional : public PoseGeneratorBase<Pose>
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

class TrajectoryDifferentiator
{
public:
    static IMUMeasurements DifferentiateTrajectory(const STLVector<InversePose>& trajectory, int deltaT_ms)
    {
        IMUMeasurements IMUs;

        if(trajectory.size() < 3)
        {
            std::cerr << "At least 3 poses are needed to differentiate a trajectory." << std::endl;
            return IMUs;
        }

        int t_ms = 0;
        IMUs.reserve(trajectory.size());
        double deltaT = 1e-3 * deltaT_ms;
        double deltaT2 = deltaT * deltaT;
        // inv_pose : local -> global, p: position in global frame, q : orientation in global frame
        size_t N = trajectory.size();
        for(size_t i = 0; i < N; i++)
        {
            Eigen::Vector3d a, omega;

            // calculate acceleration
            if(i == 0)
            {
//                a = 0.5 * (trajectory[i + 2].p - 2 * trajectory[i + 1].p + trajectory[i].p) / deltaT2;
                a = (trajectory[i + 2].p - 2 * trajectory[i + 1].p + trajectory[i].p) / deltaT2;
            }
            else if(i == 1)
            {
                a = 0.5 * (trajectory[i + 2].p - trajectory[i + 1].p - trajectory[i].p + trajectory[i - 1].p) / deltaT2;
            }
            else if(i < N - 2)
            {
                a = 0.25 * (trajectory[i + 2].p - 2 * trajectory[i].p + trajectory[i - 2].p) / deltaT2;
            }
            else if(i == N - 2)
            {
//                a = 0.5 * (trajectory[i + 1].p - 2 * trajectory[i].p + trajectory[i - 1].p) / deltaT2;
                a = 0.5 * (trajectory[i + 1].p - trajectory[i].p - trajectory[i - 1].p + trajectory[i - 2].p) / deltaT2;
            }
            else
            {
//                a = 0.5 * (trajectory[i].p - 2 * trajectory[i - 1].p + 2 * trajectory[i - 2].p) / deltaT2;
                a = (trajectory[i].p - 2 * trajectory[i - 1].p + trajectory[i - 2].p) / deltaT2;
            }

            // calculate rotational speed
            if(i < N-1)
            {
                Eigen::Quaterniond dq;
                dq.coeffs() = trajectory[i + 1].q.coeffs() - trajectory[i].q.coeffs();
                Eigen::Vector3d theta = 2 * (trajectory[i].q.inverse() * dq).vec();
                omega = theta / deltaT;
            }
            else
            {
                omega = trajectory[i].q.inverse() * (trajectory[i - 1].q * IMUs[i - 1].gyro);
            }

            Eigen::Vector3d a_local = trajectory[i].q.inverse() * a;

            IMUs.emplace_back(a_local,omega,t_ms);
            t_ms += deltaT_ms;
        }

        return IMUs;
    }
};

}  // namespace generator

#endif  // POSE_GENERATOR_FUNCTIONAL_H_