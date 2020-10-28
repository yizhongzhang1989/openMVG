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
    PoseGeneratorCircleSine(double r, double A, double fc, double fs, double T, bool storeIMU = false)
    {
        r_ = r;
        A_ = A;
        omegaCircle = 2 * PI * fc;
        omegaSine = 2 * PI * fs;
        deltaT = T;

        theta = 0;
        alpha = 0;

        timestamp = 0.0;
        storeIMU_ = storeIMU;
    }

    Pose Generate() override
    {
        Pose p;

        // generate position, t: local ->global, position of camera center, t_wc
        double x = r_ * cos(theta);
        double y = r_ * sin(theta);
        double z = A_ * sin(alpha);
        Eigen::Vector3d t(x,y,z);

        // dx = -r_ * sin(theta) * omegaCircle
        // dy = r_ * cos(theta) * omegaCircle
        // dz = A_ * cos(alpha) * omegaSine

        // dx2 = -r_ * cos(theta) * omegaCircle * omegaCircle
        // dy2 = -r_ * sin(theta) * omegaCircle * omegaCircle
        // dz2 = -A_ * sin(alpha) * omegaSine * omegaSine

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

        // generate IMU measurements
        if(storeIMU_)
        {
            double omegaCircle2 = omegaCircle * omegaCircle;
            double dx2 = -x * omegaCircle2;
            double dy2 = -y * omegaCircle2;
            double dz2 = -z * omegaSine * omegaSine;

            double sin_alpha = sin(alpha);
            double cos_alpha = cos(alpha);
            double cos2_alpha = cos_alpha * cos_alpha;
            double omegaSine2 = omegaSine * omegaSine;

            double omega_x = -omegaCircle;
            double omega_y = (sin_alpha * omegaSine2) / (1.0 + cos2_alpha * omegaSine2);
            double omega_z = 0.0;

            IMUs.emplace_back(Eigen::Vector3d(dx2,dy2,dz2),Eigen::Vector3d(omega_x,omega_y,omega_z),timestamp);
        }

        // increase position
        theta += omegaCircle * deltaT;
        alpha += omegaSine * deltaT;

        timestamp += deltaT;

        return p;
    }
    double getDeltaT() const override
    {
        return deltaT;
    }

    const IMUMeasurements& getIMUMeasurements() const
    {
        return IMUs;
    }

    bool hasIMU()
    {
        return storeIMU_;
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
    /// IMU measurements
    bool storeIMU_;
    IMUMeasurements IMUs;
    double timestamp;

};  // PoseGeneratorCircleSine

}  // namespace generator

#endif  // POSE_GENERATOR_CIRCLE_SINE_H_
