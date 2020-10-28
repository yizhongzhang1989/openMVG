#pragma once
#ifndef TYPES_H_
#define TYPES_H_

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <map>
#include <utility>
#include <vector>
#include <set>
#include <string>

namespace Eigen
{

typedef Matrix<double,6,1> Vector6d;

}  // namespace Eigen

namespace generator
{

#define PI 3.141592653589793

template<typename ValueType>
using STLVector =  std::vector<ValueType,Eigen::aligned_allocator<ValueType>>;
template<typename KeyType, typename ValueType>
using STLMap = std::map<KeyType,ValueType, std::less<KeyType>, Eigen::aligned_allocator<std::pair<const KeyType, ValueType>>>;

struct Pose
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Quaterniond q;
    Eigen::Vector3d t;
    Pose()
    {
        t.setZero();
        q.setIdentity();
    }
};

struct IMUMeasurement
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
    double timestamp;
    IMUMeasurement()
    : timestamp(0)
    {
        acc.setZero();
        gyro.setZero();
    }
    IMUMeasurement(const Eigen::Vector3d& acc_, const Eigen::Vector3d& gyro_, double t)
    : acc(acc_), gyro(gyro_), timestamp(t)
    {
        ;
    }
};
typedef STLVector<IMUMeasurement> IMUMeasurements;

struct MapPoint
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    unsigned int Id;
    Eigen::Vector3d X;
//    std::map<unsigned int, Eigen::Vector2d> obs;
    STLMap<unsigned int, Eigen::Vector2d> obs;
    MapPoint()
    :Id(-1)
    {
        obs.clear();
        X.setZero();
    }
    MapPoint(unsigned int id, const Eigen::Vector3d& X)
    : Id(id), X(X)
    {
        obs.clear();
    }
    void addObservation(unsigned int KFId, const Eigen::Vector2d& p)
    {
        obs[KFId] = p;
    }
};
//typedef std::map<unsigned int, MapPoint> MapPoints;
typedef STLMap<unsigned int, MapPoint> MapPoints;

struct KeyFrame
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    unsigned int Id;
    Pose pose;
//    std::set<unsigned int> obs;
    std::vector<unsigned int> obs;
    KeyFrame()
    :Id(-1)
    {
        obs.clear();
    }
    KeyFrame(unsigned int id, const Pose& p)
    : Id(id), pose(p)
    {
        obs.clear();
    }
    void addObservation(unsigned int PTId)
    {
//        obs.insert(PTId);
        obs.push_back(PTId);
    }
};
//typedef std::map<unsigned int, KeyFrame> KeyFrames;
typedef STLMap<unsigned int, KeyFrame> KeyFrames;

struct Simulation_Data
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    KeyFrames key_frames;
    MapPoints map_points;
    Simulation_Data()
    {
        key_frames.clear();
        map_points.clear();
    }
};

}  // namespace generator

#endif  // TYPES_H_