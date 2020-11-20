//
// Created by root on 10/16/20.
//

#ifndef OPENMVG_VISFM_CERES_PARAM_HPP
#define OPENMVG_VISFM_CERES_PARAM_HPP

//#include <eigen3/Eigen/Dense>
#include <Eigen/Dense>
#include <ceres/ceres.h>
//#include "third_party/ceres-solver/include/ceres/ceres.h"
#include "Utility.hpp"



class PoseLocalParameterization : public ceres::LocalParameterization
{
    virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const;
    virtual bool ComputeJacobian(const double *x, double *jacobian) const;
    virtual int GlobalSize() const { return 7; };
    virtual int LocalSize() const { return 6; };
};

class PoseQuaternLocalParameterization : public ceres::LocalParameterization
{
    virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const;
    virtual bool ComputeJacobian(const double *x, double *jacobian) const;
    virtual int GlobalSize() const { return 4; };
    virtual int LocalSize() const { return 3; };
};

#endif //OPENMVG_VISFM_CERES_PARAM_HPP
