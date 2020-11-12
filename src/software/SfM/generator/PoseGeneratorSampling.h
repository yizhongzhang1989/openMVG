#pragma once
#ifndef POSE_GENERATOR_SAMPLING_H_
#define POSE_GENERATOR_SAMPLING_H_

#include <utility>

#include "PoseGeneratorBase.h"
#include "types.h"
#include "TrajectorySampler.h"
#include "TrajectoryDifferentiator.h"

namespace generator
{

class PoseGeneratorSampling : public PoseGeneratorBase<Pose, Eigen::aligned_allocator<Pose>>
{
public:
    PoseGeneratorSampling(std::string  sFileName, double total_duration_s, int deltaT, int deltaT_IMU, bool storeIMU = false, LookDirection direction = FORWARD)
    : mFileName(std::move(sFileName)), total_duration(total_duration_s), direction(direction), storeIMU_(storeIMU), deltaT(deltaT), deltaT_IMU(deltaT_IMU)
    {
        double T_sample_pose = deltaT * 1e-3;
        mInvPoses = TrajectorySampler::SampleTrajectory(mFileName,total_duration,T_sample_pose);
        mPoses.clear();
        mPoses.reserve(mInvPoses.size());
        for(const InversePose & inv_pose : mInvPoses)
        {
            Pose pose;
            pose.q = inv_pose.q.inverse();
            pose.t = -(pose.q * inv_pose.p);
            mPoses.push_back(pose);
        }

        if(storeIMU_)
        {
            double T_sample_imu = deltaT_IMU * 1e-3;
            mInvPosesIMU = TrajectorySampler::SampleTrajectory(mFileName,total_duration,T_sample_imu);
            TotalIMUs = TrajectoryDifferentiator::DifferentiateTrajectory(mInvPosesIMU,deltaT_IMU);
        }
        else
        {
            TotalIMUs.clear();
            mInvPosesIMU.clear();
        }
        IMUs.clear();
    }

    Pose Generate() override
    {
        throw std::runtime_error("[PoseGeneratorSampling] Method Generate() is disabled.");
        return Pose();
    }

    Poses Generate(int num_poses) override
    {
        if(num_poses == -1)
        {
            if(storeIMU_)
            {
                IMUs = TotalIMUs;
            }

            return mPoses;
        }

        Poses poses;
        int t_cam = -deltaT;
        for(int i = 0; i < num_poses && i < mPoses.size(); i++)
        {
            poses.push_back(mPoses[i]);
            t_cam += deltaT;
        }

        if(storeIMU_)
        {
            IMUs.clear();
            int t_imu = 0;
            int idx_imu = 0;
            while(t_imu <= t_cam)
            {
                IMUs.push_back(TotalIMUs[idx_imu]);
                idx_imu++;
                t_imu += deltaT_IMU;
            }
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

    const STLVector<InversePose>& getInversePoses() const
    {
        return mInvPoses;
    }

    const STLVector<InversePose>& getInversePosesIMU() const
    {
        return mInvPosesIMU;
    }

private:
    LookDirection direction;
    // sampling period in ms
    int deltaT;
    int deltaT_IMU;
    // IMU measurements
    IMUMeasurements TotalIMUs;
    IMUMeasurements IMUs;
    const bool storeIMU_;
    // for trajectory sampling
    const std::string mFileName;
    const double total_duration;
    STLVector<Pose> mPoses;
    STLVector<InversePose> mInvPoses;
    STLVector<InversePose> mInvPosesIMU;
};

}  // namespace generator

#endif  // POSE_GENERATOR_SAMPLING_H_