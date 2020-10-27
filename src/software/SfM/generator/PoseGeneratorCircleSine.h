#pragma once
#ifndef POSE_GENERATOR_CIRCLE_SINE_H_
#define POSE_GENERATOR_CIRCLE_SINE_H_

#include "PoseGeneratorBase.h"
#include "types.h"

namespace generator
{

class PoseGeneratorCircleSine : public PoseGeneratorBase<Pose>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// r: radius of circle;  A: amptitude of sine;
    /// fc: rotation frequency (r/s);  fs: frequency of sine;
    /// T: sampling period.
    PoseGeneratorCircleSine(double r, double A, double fc, double fs, double T)
    {
        r_ = r;
        A_ = A;
        omegaCircle = 2 * PI * fc;
        omegaSine = 2 * PI * fs;
        deltaT = T;

        theta = 0;
        alpha = 0;
    }

    Pose Generate() override
    {
        Pose p;

        // generate position, t: local ->global, position of camera center, t_wc
        double x = r_ * cos(theta);
        double y = r_ * sin(theta);
        double z = A_ * sin(alpha);
        Eigen::Vector3d t(x,y,z);

        // generate orientation
        Eigen::Vector3d rz(-y*omegaCircle,x*omegaCircle,A_*cos(alpha)*omegaSine);
        Eigen::Vector3d ry(-x,-y,0);
        Eigen::Vector3d rx = ry.cross(rz);

        rx = rx / rx.norm();
        ry = ry / ry.norm();
        rz = rz / rz.norm();

        // R: local -> global, R_wc
        Eigen::Matrix3d R;
        R.col(0) = rx;
        R.col(1) = ry;
        R.col(2) = rz;

        // p.q: global -> local, q_cw
        p.q = Eigen::Quaterniond(R.transpose());
        // p.t: global -> local, t_cw
        p.t = - R.transpose() * t;

        // increase position
        theta += omegaCircle * deltaT;
        alpha += omegaSine * deltaT;

        return p;
    }
    double getDeltaT() const override
    {
        return deltaT;
    }

private:
    double r_;
    double A_;
    double omegaCircle;
    double omegaSine;
    double deltaT;
    /// angle of rotation (rad)
    double theta;
    /// angle of sine (rad)
    double alpha;

};  // PoseGeneratorCircleSine

}  // namespace generator

#endif  // POSE_GENERATOR_CIRCLE_SINE_H_
