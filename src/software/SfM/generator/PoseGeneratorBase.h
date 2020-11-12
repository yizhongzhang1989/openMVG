#pragma once
#ifndef POSE_GENERATOR_BASE_H_
#define POSE_GENERATOR_BASE_H_

#include <vector>

namespace generator
{

#define GRAVITY 0.0

template<class PoseT, class Allocator = std::allocator<PoseT>>
class PoseGeneratorBase
{
public:
    enum LookDirection
    {
        FORWARD,
        LEFTWARD
    };
    typedef PoseT pose_type;
    typedef std::vector<PoseT,Allocator> Poses;
public:
    // generate single pose
    virtual pose_type Generate() = 0;
    // generate a number of poses
    virtual Poses Generate(int num_poses) = 0;
    // get sampling period in ms
    virtual int getDeltaT() const = 0;
};

}  // namespace generator

#endif  // POSE_GENERATOR_BASE_H_
