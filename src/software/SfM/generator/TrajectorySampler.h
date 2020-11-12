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
    static STLVector<InversePose> SampleTrajectory(const std::string& objFile, double total_duration, double T_sample, STLVector<Eigen::Vector3d>* pOriginalVerticesOut = nullptr)
    {
        // read squares
        STLVector<Eigen::Vector3d> vertices;
        ObjFileReader::ReadVerticesOnlyFromObj(objFile.c_str(),vertices);

        std::cout<<"num vertices: "<<vertices.size()<<std::endl;

        if(pOriginalVerticesOut)
            *pOriginalVerticesOut = vertices;

        int N = vertices.size() / 2;
        STLVector<Eigen::Vector3d> up_vertices, down_vertices;

        for(int i = 0; i < N; i++)
        {
            down_vertices.push_back(vertices[2 * i]);
            up_vertices.push_back(vertices[2 * i + 1]);
        }

        // interpolate grid
        STLVector<Eigen::Vector3d> up_vertices_interpolated, down_vertices_interpolated;
        InterpolateGrid(up_vertices,down_vertices,total_duration,T_sample,up_vertices_interpolated,down_vertices_interpolated);

        // generate poses
        STLVector<InversePose> inv_poses;
        GeneratePosesFromGrid(up_vertices_interpolated,down_vertices_interpolated,inv_poses);

        return inv_poses;
    }
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
        GeneratePosesFromGrid(up_vertices,down_vertices,inv_poses);

        return inv_poses;
    }
private:
    static void GeneratePosesFromGrid(const STLVector<Eigen::Vector3d>& up_vertices, const STLVector<Eigen::Vector3d>& down_vertices, STLVector<InversePose>& inv_poses)
    {
        assert(up_vertices.size() == down_vertices.size());

        inv_poses.clear();
        int N = up_vertices.size();
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
    }
    static void InterpolateGrid(const STLVector<Eigen::Vector3d>& up_vertices_in, const STLVector<Eigen::Vector3d>& down_vertices_in,
                         const double total_duration, const double T_sample,
                         STLVector<Eigen::Vector3d>& up_vertices_out, STLVector<Eigen::Vector3d>& down_vertices_out)
    {
        assert(up_vertices_in.size() == down_vertices_in.size());
        assert(total_duration > 0 && T_sample > 0);

        int N_vertex = up_vertices_in.size();
        int N_grid = N_vertex - 1;
        double T_grid = total_duration / N_grid;

        up_vertices_out.clear();
        down_vertices_out.clear();

        double current_time = 0.0;
        while(current_time < total_duration)
        {
            double pos = current_time / T_grid;
            int grid_i = int(pos);
            double remain = pos - grid_i;

            // linear interpolation
            const Eigen::Vector3d& A = up_vertices_in[grid_i];
            const Eigen::Vector3d& B = up_vertices_in[grid_i + 1];
            const Eigen::Vector3d& C = down_vertices_in[grid_i];
            const Eigen::Vector3d& D = down_vertices_in[grid_i + 1];

            Eigen::Vector3d up_vertex = A + remain* (B - A);
            Eigen::Vector3d down_vertex = C + remain * (D - C);

            up_vertices_out.push_back(up_vertex);
            down_vertices_out.push_back(down_vertex);

            current_time += T_sample;
        }
    }
};





}  // namespace generator

#endif  // TRAJECTORY_SAMPLER_H_