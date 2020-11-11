#pragma once
#ifndef TRAJECTORY_SAMPLER_H_
#define TRAJECTORY_SAMPLER_H_

#include <iostream>
#include "OBJFileReader.h"
#include "types.h"

namespace generator
{

class TrajectorySampler
{
public:
    static STLVector<InversePose> SampleTrajectory(const std::string& objFile, STLVector<Eigen::Vector3d>* pVerticesOut = nullptr)
    {
        // read squares
        STLVector<Eigen::Vector3d> vertices;
        ObjFileReader::ReadVerticesOnlyFromObj(objFile.c_str(),vertices);

        std::cout<<"num vertices: "<<vertices.size()<<std::endl;

        if(pVerticesOut)
            *pVerticesOut = vertices;

        int N = vertices.size() / 2;
        STLVector<Eigen::Vector3d> up_vertices, down_vertices;

        for(int i = 0; i < N; i++)
        {
            down_vertices.push_back(vertices[2 * i]);
            up_vertices.push_back(vertices[2 * i + 1]);
        }

        // generate poses
        STLVector<InversePose> inv_poses;
        for(int i = 0; i < N - 1; i++)
        {
            InversePose inv_pose;
            inv_pose.p = down_vertices[i];

            Eigen::Vector3d OA = down_vertices[i + 1] - down_vertices[i];
            Eigen::Vector3d OC = up_vertices[i] - down_vertices[i];

            Eigen::Vector3d rz = OA;
            Eigen::Vector3d rx = OA.cross(OC);
            Eigen::Vector3d ry = rz.cross(rx);

            rx.normalize();
            ry.normalize();
            rz.normalize();

            Eigen::Matrix3d R;
            R.col(0) = rx;
            R.col(1) = ry;
            R.col(2) = rz;

            inv_pose.q = Eigen::Quaterniond(R);

            inv_poses.push_back(inv_pose);
        }

        return inv_poses;
    }
};





}  // namespace generator

#endif  // TRAJECTORY_SAMPLER_H_