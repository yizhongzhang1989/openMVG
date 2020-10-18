//
// Created by root on 10/16/20.
//

#ifndef OPENMVG_VISFM_CERES_PARAM_H
#define OPENMVG_VISFM_CERES_PARAM_H

#include <eigen3/Eigen/Dense>
#include <ceres/ceres.h>
#include "Utility.hpp"



class PoseLocalParameterization : public ceres::LocalParameterization
{
    virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const;
    virtual bool ComputeJacobian(const double *x, double *jacobian) const;
    virtual int GlobalSize() const { return 7; };
    virtual int LocalSize() const { return 6; };
};

#endif //OPENMVG_VISFM_CERES_PARAM_H
