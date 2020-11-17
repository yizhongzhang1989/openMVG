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
        std::cout << "total duration = " << total_duration << " , T_sample = " << T_sample << std::endl;
        std::cout << up_vertices_interpolated.size() << " " << down_vertices_interpolated.size() << std::endl;

        // generate poses
        STLVector<InversePose> inv_poses;
        GeneratePosesFromGrid(up_vertices_interpolated,down_vertices_interpolated,inv_poses);
        CorrectQuaternion(inv_poses);

        std::cout << "num_poses = " << inv_poses.size() << std::endl;
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

            //  =======================================
            //  CatmullRom spline interpolation
			Eigen::Vector3d up_v0 = grid_i > 0 ? up_vertices_in[grid_i - 1] : (up_vertices_in[grid_i] - up_vertices_in[grid_i + 1]) * 3 + up_vertices_in[grid_i + 2]; /*up_vertices_in[grid_i] * 2 - up_vertices_in[grid_i + 1];*/
			Eigen::Vector3d up_v1 = up_vertices_in[grid_i];
			Eigen::Vector3d up_v2 = up_vertices_in[grid_i + 1];
			Eigen::Vector3d up_v3 = grid_i < N_grid - 1 ? up_vertices_in[grid_i + 2] : up_vertices_in[grid_i + 1] * 2 - up_vertices_in[grid_i];
			Eigen::Vector3d up_vertex = InterpCatmullRom(up_v0, up_v1, up_v2, up_v3, remain);

			Eigen::Vector3d down_v0 = grid_i > 0 ? down_vertices_in[grid_i - 1] : (down_vertices_in[grid_i] - down_vertices_in[grid_i + 1]) * 3 + down_vertices_in[grid_i + 2];/*down_vertices_in[grid_i] * 2 - down_vertices_in[grid_i + 1];*/
            Eigen::Vector3d down_v1 = down_vertices_in[grid_i];
            Eigen::Vector3d down_v2 = down_vertices_in[grid_i + 1];
            Eigen::Vector3d down_v3 = grid_i < N_grid - 1 ? down_vertices_in[grid_i + 2] : down_vertices_in[grid_i + 1] * 2 - down_vertices_in[grid_i];
            Eigen::Vector3d down_vertex = InterpCatmullRom(down_v0, down_v1, down_v2, down_v3, remain);


            //  =======================================
            //// linear interpolation
            //const Eigen::Vector3d& A = up_vertices_in[grid_i];
            //const Eigen::Vector3d& B = up_vertices_in[grid_i + 1];
            //const Eigen::Vector3d& C = down_vertices_in[grid_i];
            //const Eigen::Vector3d& D = down_vertices_in[grid_i + 1];

            //Eigen::Vector3d up_vertex = A + remain* (B - A);
            //Eigen::Vector3d down_vertex = C + remain * (D - C);


            //  =======================================
            up_vertices_out.push_back(up_vertex);
            down_vertices_out.push_back(down_vertex);

            current_time += T_sample;
        }

		//static int flag = 0;
		//char filename[1024];
		//sprintf(filename, "E:/OpenMVG_IMU_synthetic_data/circle/tmp_%d.txt", flag);
		//std::ofstream file(filename);
		//for (int i = 0; i < down_vertices_out.size(); i++) {
		//	Eigen::Vector3d r = up_vertices_out[i];
		//	file << r[0] << ',' << r[1] << ',' << r[2] << std::endl;
		//}
		//file.close();
		//flag++;
	}

	template<typename T>
	static T InterpCatmullRom(T val0, T val1, T val2, T val3, double t) {
		double t2 = t * t;

		return
			(val0 * (-0.5) + val1 * 1.5 + val2 * (-1.5) + val3 * (0.5)) * t2 * t +
			(val0 + val1 * (-2.5) + val2 * 2 + val3 * (-0.5)) * t2 +
			(val0 * (-0.5) + val2 * (0.5)) * t + val1;
	}

    static void CorrectQuaternion(STLVector<InversePose>& trajectory)
    {
        if (trajectory.empty())
            return;
        Eigen::Quaterniond last_q = trajectory[0].q;
        for (InversePose& inv_pose : trajectory)
        {
            Eigen::Quaterniond& q = inv_pose.q;
            Eigen::Quaterniond tmp_q;
            tmp_q.coeffs() = -q.coeffs();

            double norm1 = (q.coeffs() - last_q.coeffs()).norm();
            double norm2 = (tmp_q.coeffs() - last_q.coeffs()).norm();
            if (norm1 > norm2)
            {
                q = tmp_q;
            }

            last_q = q;
        }
    }

};





}  // namespace generator

#endif  // TRAJECTORY_SAMPLER_H_