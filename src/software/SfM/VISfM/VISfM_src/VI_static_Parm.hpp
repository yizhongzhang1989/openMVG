//
// Created by xin on 2020/10/18.
//

#ifndef OPENMVG_VI_STATIC_PARM_HPP
#define OPENMVG_VI_STATIC_PARM_HPP

//#include <eigen3/Eigen/Dense>
#include <Eigen/Dense>

class VIstaticParm
{
public:
    static Eigen::Vector3d G_;
    static double gyr_n;
    static double gyr_w;
    static double acc_n;
    static double acc_w;
};

#endif //OPENMVG_VI_STATIC_PARM_HPP
