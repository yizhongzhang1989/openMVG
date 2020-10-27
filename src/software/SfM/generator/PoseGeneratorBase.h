#pragma once
#ifndef POSE_GENERATOR_BASE_H_
#define POSE_GENERATOR_BASE_H_

namespace generator
{

template<typename PoseT>
class PoseGeneratorBase
{
public:
    typedef PoseT pose_type;
    virtual pose_type Generate() = 0;
    virtual double getDeltaT() const = 0;
};

}  // namespace generator

#endif  // POSE_GENERATOR_BASE_H_
