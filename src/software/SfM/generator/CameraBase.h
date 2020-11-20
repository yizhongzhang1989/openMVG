#pragma once
#ifndef CAMERA_BASE_H_
#define CAMERA_BASE_H_

#include <string>
#include <vector>

namespace generator
{

template<typename Point3D, typename Point2D>
class CameraBase
{
public:
    virtual Point2D Project(const Point3D&) const = 0;
    virtual bool isInView(const Point2D&) const = 0;
    virtual std::string getModelName() const = 0;
    virtual int getWidth() const = 0;
    virtual int getHeight() const = 0;
    virtual void getParams(std::vector<double>&) const = 0;
};  // CameraBase

}  // namespace generator

#endif  // CAMERA_BASE_H_