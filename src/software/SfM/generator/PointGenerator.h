#pragma once
#ifndef POINT_GENERATOR_H_
#define POINT_GENERATOR_H_

#include <Eigen/Core>
#include "PointGeneratorBase.h"
#include "types.h"
#include "SurfaceSampler.h"

namespace generator
{

class PointGenerator : public PointGeneratorBase<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PointGenerator(double minx, double maxx, double miny, double maxy, double minz, double maxz)
    : PointGeneratorBase(minx, maxx, miny, maxy, minz, maxz)
    {
        ;
    }
    explicit PointGenerator(const std::string& sFileName)
    : PointGeneratorBase(sFileName)
    {
        ;
    }
private:
    bool GenerateSampling(const std::string& sFileName, int num_points, Points& points) override
    {
        points = TriangleSampler::SampleObjFile(sFileName.c_str(),num_points);
        return true;
    }
};

}  // namespace generator

#endif  // POINT_GENERATOR_H_