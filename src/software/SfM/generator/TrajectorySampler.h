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
    //  cubic spline interpolation with natrual boundary condition
    template <class T>
    class InterpCubicSplineNatural {
    public:
        InterpCubicSplineNatural(const std::vector<double>& x, const STLVector<T>& y) {
            //	validate data
            if (x.size() < 4) {
                std::cout << "error: InterpCubicSplineNatural, too few points" << std::endl;
                return;
            }
            if (x.size() != y.size()) {
                std::cout << "error: InterpCubicSplineNatural, x.size() != y.size()" << std::endl;
                return;
            }
            for (int i = 1; i < x.size(); i++) {
                if (x[i] <= x[i - 1]) {
                    std::cout << "error: InterpCubicSplineNatural, x must be increasing order" << std::endl;
                    return;
                }
            }

            this->x = x;
            this->y = y;

            PartialDirivitive();
        }

        T Interp(double x) {
            if (x <= this->x[0])
                return this->y[0];
            if (x >= this->x.back())
                return this->y.back();

            auto iter = std::upper_bound(this->x.begin(), this->x.end(), x);
            int idx = iter - this->x.begin() - 1;
            if (idx < 0 || idx >= a.size()) {
                std::cout << "error: InterpCubicSplineNatural::Interp, idx not correct: " << idx << std::endl;
                return this->y[0];
            }

            double coef = x - this->x[idx];
            double coef2 = coef * coef;
            return a[idx] + b[idx] * coef + c[idx] * coef2 + d[idx] * coef2 * coef;
        }

    protected:
        std::vector<double>		x;
        STLVector<T>			y;

        std::vector<T>			a, b, c, d;		//	a + b (x-xi) + c (x-xi)^2 + d (x-xi)^3

        void PartialDirivitive() {
            //	the algorithm can be found here: https://blog.csdn.net/limingmcu/article/details/91492214
            std::vector<double> h(x.size() - 1, 0.0);
            for (int i = 0; i < h.size(); i++)
                h[i] = x[i + 1] - x[i];

            //	setup the tri-digonal matrix
            std::vector<double> a(x.size(), 0.0), b(x.size(), 0.0), c(x.size(), 0.0);
            b[0] = 1.0;
            b[x.size() - 1] = 1.0;
            for (int i = 1; i < x.size() - 1; i++) {
                a[i] = h[i - 1];
                b[i] = 2 * (h[i - 1] + h[i]);
                c[i] = h[i];
            }

            //	setup matrix right hand side
            T zero = Zero();
            std::vector<T> rhs(x.size(), zero);
            for (int i = 1; i < x.size() - 1; i++) {
                rhs[i] = ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]) * 6.0;
            }

            //	solve the system using Thomas method: https://wenku.baidu.com/view/5503af3589eb172dec63b719.html
            std::vector<double> l(x.size(), 0.0), u(x.size(), 0.0);
            u[0] = b[0];
            for (int i = 1; i < x.size(); i++) {
                l[i] = a[i] / u[i - 1];
                u[i] = b[i] - c[i - 1] * l[i];
            }

            std::vector<T> f(x.size(), zero);
            f[0] = rhs[0];
            for (int i = 1; i < x.size(); i++) {
                f[i] = rhs[i] - l[i] * f[i - 1];
            }

            std::vector<T> m(x.size(), zero);
            m[x.size() - 1] = f[x.size() - 1] / u[x.size() - 1];
            for (int i = x.size() - 2; i >= 0; i--) {
                m[i] = (f[i] - c[i] * m[i + 1]) / u[i];
            }

            //	create third order coef for each segment
            this->a.resize(x.size() - 1, zero);
            this->b.resize(x.size() - 1, zero);
            this->c.resize(x.size() - 1, zero);
            this->d.resize(x.size() - 1, zero);
            for (int i = 0; i < x.size() - 1; i++) {
                this->a[i] = y[i];
                this->b[i] = (y[i + 1] - y[i]) / h[i] - h[i] / 2.0 * m[i] - h[i] / 6.0 * (m[i + 1] - m[i]);
                this->c[i] = m[i] / 2.0;
                this->d[i] = (m[i + 1] - m[i]) / (6.0 * h[i]);
            }
        }

        T Zero() {
            return T();
        }
    };


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

        std::vector<double> coef;
        for (int i = 0; i < N_vertex; i++)
			coef.push_back(i);

		InterpCubicSplineNatural<Eigen::Vector3d> cubic_spline_interp_up(coef, up_vertices_in);
        InterpCubicSplineNatural<Eigen::Vector3d> cubic_spline_interp_down(coef, down_vertices_in);

        double current_time = 0.0;
        while (current_time < total_duration)
        {
            double pos = current_time / T_grid;

            Eigen::Vector3d up_vertex = cubic_spline_interp_up.Interp(pos);
            Eigen::Vector3d down_vertex = cubic_spline_interp_down.Interp(pos);

            //  =======================================
            up_vertices_out.push_back(up_vertex);
            down_vertices_out.push_back(down_vertex);

            current_time += T_sample;
        }

        //static int flag = 0;
        //char filename[1024];
        //sprintf(filename, "E:/OpenMVG_IMU_synthetic_data/circle_y/tmp_%d.txt", flag);
        //std::ofstream file(filename);
        //for (int i = 0; i < down_vertices_out.size(); i++) {
        //	Eigen::Vector3d r = up_vertices_out[i];
        //	file << r[0] << ',' << r[1] << ',' << r[2] << std::endl;
        //}
        //file.close();
        //flag++;
    }

    static void InterpolateGridFirstOrder(const STLVector<Eigen::Vector3d>& up_vertices_in, const STLVector<Eigen::Vector3d>& down_vertices_in,
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