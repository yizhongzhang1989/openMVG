#pragma once
#ifndef POINT_GENERATOR_H_
#define POINT_GENERATOR_H_

#include <random>
#include <chrono>
#include <Eigen/Core>

namespace generator
{

class PointGenerator
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    PointGenerator(double minx, double maxx, double miny, double maxy, double minz, double maxz)
    :e(std::chrono::system_clock::now().time_since_epoch().count())
    {
        nx = std::uniform_real_distribution<double>(minx,maxx);
        ny = std::uniform_real_distribution<double>(miny,maxy);
        nz = std::uniform_real_distribution<double>(minz,maxz);
    }
    Eigen::Vector3d Generate()
    {
        double x = nx(e);
        double y = ny(e);
        double z = nz(e);

        return {x,y,z};
    }
private:
    std::default_random_engine e;
    std::uniform_real_distribution<double> nx,ny,nz;
};  // class PointGenerator

}  // namespace generator

#endif  // POINT_GENERATOR_H_
