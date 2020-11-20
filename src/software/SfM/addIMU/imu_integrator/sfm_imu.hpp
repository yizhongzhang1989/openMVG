// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_IMU_HPP
#define OPENMVG_SFM_SFM_IMU_HPP

#include <vector>
#include <string>
#include <algorithm>

#include <Eigen/Core>

#include "openMVG/types.hpp"

namespace openMVG {
namespace sfm {

template<typename Data_T>
struct SfM_IMU_Base
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using accelerometer_t = Eigen::Matrix<Data_T,3,1>;
    using gyroscope_t = Eigen::Matrix<Data_T,3,1>;
    virtual void ReadFromCSV(const std::string& filename) = 0;
protected:
    std::vector<accelerometer_t,Eigen::aligned_allocator<accelerometer_t>> vec_acc;
    std::vector<gyroscope_t,Eigen::aligned_allocator<gyroscope_t>> vec_gyro;
    std::vector<double> timestamps;

    bool checkTimestamp(double timestamp, double eps)
    {
        if(timestamps.empty())
            return false;
        if(timestamps[0]-timestamp>eps || timestamp-timestamps.back()>eps)
            return false;
        return true;
    }
public:
    accelerometer_t getAcceleration(IndexT idx) const
    {
        assert(idx>=0 && idx<vec_acc.size());
        return vec_acc[idx];
    }
    gyroscope_t getOmega(IndexT idx) const
    {
        assert(idx>=0 && idx<vec_gyro.size());
        return vec_gyro[idx];
    }
    accelerometer_t getAcceleration(double timestamp) const
    {
        assert(checkTimestamp(timestamp, 0.1));
        auto it = std::lower_bound(timestamps.begin(),timestamps.end(),timestamp);

        if(it == timestamps.begin())
            return vec_acc[0];
        else if(it == timestamps.end())
            return vec_acc.back();

        IndexT idx = it-timestamps.begin();
        if(fabs(timestamps[idx-1]-timestamp)<fabs(timestamps[idx]-timestamp))
            return vec_acc[idx-1];
        return vec_acc[idx];
    }
    gyroscope_t getOmega(double timestamp) const
    {
        assert(checkTimestamp(timestamp, 0.1));
        auto it = std::lower_bound(timestamps.begin(),timestamps.end(),timestamp);

        if(it == timestamps.begin())
            return vec_gyro[0];
        else if(it == timestamps.end())
            return vec_gyro.back();

        IndexT idx = it-timestamps.begin();
        if(fabs(timestamps[idx-1]-timestamp)<fabs(timestamps[idx]-timestamp))
            return vec_gyro[idx-1];
        return vec_gyro[idx];
    }
    size_t Size() const
    {
        return vec_acc.size();
    }
};

struct SfM_IMU : public SfM_IMU_Base<double>
{
    virtual void ReadFromCSV(const std::string& filename);
    SfM_IMU()
    {
        vec_acc.clear();
        vec_gyro.clear();
        timestamps.clear();
    }
};

}  // namespace sfm
}  // namespace openMVG

#endif  // OPENMVG_SFM_SFM_IMU_HPP