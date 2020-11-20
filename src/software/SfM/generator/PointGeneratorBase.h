#pragma once
#ifndef POINT_GENERATOR_BASE_H_
#define POINT_GENERATOR_BASE_H_

#include <random>
#include <chrono>
#include <utility>
#include <vector>

namespace generator
{

template<class PointType, class Allocator = std::allocator<PointType>>
class PointGeneratorBase
{
public:
    enum PointGenerationMode
    {
        RANDOM_POINT,
        SAMPLING_SURFACE
    };
    typedef PointType Point;
    typedef std::vector<PointType,Allocator> Points;

    PointGeneratorBase(double minx, double maxx, double miny, double maxy, double minz, double maxz)
    :e(std::chrono::system_clock::now().time_since_epoch().count()), mPointMode(RANDOM_POINT)
    {
        nx = std::uniform_real_distribution<double>(minx,maxx);
        ny = std::uniform_real_distribution<double>(miny,maxy);
        nz = std::uniform_real_distribution<double>(minz,maxz);
    }
    explicit PointGeneratorBase(std::string sFileName)
    : mFileName(std::move(sFileName)), mPointMode(SAMPLING_SURFACE)
    {
        ;
    }

    Points Generate(int num_points)
    {
        Points points;

        switch(mPointMode)
        {
            case RANDOM_POINT:
            {
                for(int i = 0; i < num_points; i++)
                {
                    points.push_back(GenerateRandom());
                }
            }break;
            case SAMPLING_SURFACE:
            {
                if(!GenerateSampling(mFileName,num_points,points))
                {
                    throw std::runtime_error("[PointGenerator] sampling surface failed.");
                }
            }break;
        }

        return points;
    }

private:
    // generation mode
    const PointGenerationMode mPointMode;

    // for random points
    std::default_random_engine e;
    std::uniform_real_distribution<double> nx,ny,nz;

    // for sampling surface
    const std::string mFileName;

    Point GenerateRandom()
    {
        if(mPointMode != RANDOM_POINT)
            throw std::runtime_error("[PointGenerator] random generator is not setup.");

        double x = nx(e);
        double y = ny(e);
        double z = nz(e);

        return {x,y,z};
    }

    virtual bool GenerateSampling(const std::string& sFileName, int num_points, Points& points)
    {
        return false;
    }

};  // class PointGenerator

}  // namespace generator

#endif  // POINT_GENERATOR_BASE_H_
