#pragma once
#ifndef TRAJECTORY_DIFFERENTIATOR_H_
#define TRAJECTORY_DIFFERENTIATOR_H_

#include "types.h"

namespace generator
{

class TrajectoryDifferentiator
{
public:
    static IMUMeasurements DifferentiateTrajectory(const STLVector<InversePose>& trajectory_original, int deltaT_ms, bool correct_quaternion = true)
    {
        IMUMeasurements IMUs;

        if(trajectory_original.size() < 3)
        {
            std::cerr << "At least 3 poses are needed to differentiate a trajectory." << std::endl;
            return IMUs;
        }

        const STLVector<InversePose>* trajectory_ptr = nullptr;
        STLVector<InversePose> trajectory_buff;
        if(correct_quaternion)
        {
            trajectory_buff = trajectory_original;
            CorrectQuaternion(trajectory_buff);
            trajectory_ptr = &trajectory_buff;
        }
        else
        {
            trajectory_ptr = &trajectory_original;
        }
        const STLVector<InversePose>& trajectory = *trajectory_ptr;

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
private:
    static void CorrectQuaternion(STLVector<InversePose>& trajectory)
    {
        if (trajectory.empty())
            return;
        Eigen::Quaterniond last_q = trajectory[0].q;
        for (InversePose& inv_pose : trajectory)
        {
            Eigen::Quaterniond& q = inv_pose.q;
            Eigen::Quaterniond tmp_q;
            tmp_q.coeffs() = -q.coeffs();

            double norm1 = (q.coeffs() - last_q.coeffs()).norm();
            double norm2 = (tmp_q.coeffs() - last_q.coeffs()).norm();
            if (norm1 > norm2)
            {
                q = tmp_q;
            }

            last_q = q;
        }
    }
};



}


#endif