#pragma once
#ifndef POSE_GENERATOR_BASE_H_
#define POSE_GENERATOR_BASE_H_

namespace generator
{

#define GRAVITY 0.0

template<typename PoseT>
class PoseGeneratorBase
{
public:
    enum LookDirection
    {
        FORWARD,
        LEFTWARD
    };
    typedef PoseT pose_type;
    virtual pose_type Generate() = 0;
    // get sampling period in ms
    virtual int getDeltaT() const = 0;
};

}  // namespace generator

#endif  // POSE_GENERATOR_BASE_H_
