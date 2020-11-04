#pragma once
#ifndef CAMERA_PINHOLE_H_
#define CAMERA_PINHOLE_H_

#include "CameraBase.h"
#include <Eigen/Core>

namespace generator
{

class CameraPinhole : public CameraBase<Eigen::Vector3d,Eigen::Vector2d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    CameraPinhole(double fx, double fy, double cx, double cy, int width, int height)
    :fx_(fx), fy_(fy), cx_(cx), cy_(cy), width_(width), height_(height)
    {
        ;
    }
    Eigen::Vector2d Project(const Eigen::Vector3d& X) const override
    {
        double xp = fx_*X[0]/X[2]+cx_;
        double yp = fy_*X[1]/X[2]+cy_;
        return {xp, yp};
    }
    bool isInView(const Eigen::Vector2d& p) const override
    {
        return (p[0]>=0 && p[0]<width_ && p[1]>=0 && p[1]<height_);
    }
    std::string getModelName() const override
    {
        return "PINHOLE";
    }
    int getWidth() const override
    {
        return width_;
    }
    int getHeight() const override
    {
        return height_;
    }
    void getParams(std::vector<double>& params) const override
    {
        params.clear();
        params.push_back(fx_);
        params.push_back(fy_);
        params.push_back(cx_);
        params.push_back(cy_);
    }
private:
    double fx_,fy_,cx_,cy_;
    int width_,height_;
};

}  // namespace generator

#endif  // CAMERA_PINHOLE_H_