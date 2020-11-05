#pragma once
#ifndef SURFACE_SAMPLER_H_
#define SURFACE_SAMPLER_H_

#include <random>
#include <chrono>
#include "types.h"
#include "OBJFileReader.h"

namespace generator
{

class TriangleSampler
{
public:
    static STLVector<Eigen::Vector3d> SampleTriangle(const Triangle& triangle, int num_points)
    {
        std::default_random_engine e(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> n(0.0,1.0);
        STLVector<Eigen::Vector3d> samples;
        for(int i=0;i<num_points;i++)
        {
            double a = n(e);
            double b = n(e);
            a = sqrt(a);
            Eigen::Vector3d p = (1 - a) * triangle.v1 + a * (1 - b) * triangle.v2 + a * b * triangle.v3;
            samples.push_back(p);
        }
        return samples;
    }
    static double CalculateTriangleArea(const Triangle& triangle)
    {
        double a = (triangle.v1 - triangle.v2).norm();
        double b = (triangle.v1 - triangle.v3).norm();
        double c = (triangle.v2 - triangle.v3).norm();
        double p = 0.5 * (a + b + c);
        double area = p * (p - a) * (p - b) * (p - c);
        area = sqrt(area);
        return area;
    }
    static STLVector<Eigen::Vector3d> SampleMultipleTriangles(const STLVector<Triangle>& triangles, int num_points)
    {
        std::vector<double> areas;
        STLVector<Eigen::Vector3d> point_cloud;
        double sum_area = 0.0;
        for(const Triangle& triangle : triangles)
        {
            double area = CalculateTriangleArea(triangle);
            areas.push_back(area);
            sum_area += area;
        }
        for(size_t i = 0; i < triangles.size(); i++)
        {
            int num_points_i = std::ceil(areas[i] * num_points / sum_area);
            STLVector<Eigen::Vector3d> point_cloud_i = SampleTriangle(triangles[i], num_points_i);
            point_cloud.insert(point_cloud.end(),point_cloud_i.begin(),point_cloud_i.end());
        }
        return point_cloud;
    }
    static STLVector<Eigen::Vector3d> SampleMultipleTrianglesWithPos(const STLVector<Triangle>& triangles, int num_points, const Eigen::Vector3d& position, bool visible_only = false)
    {
        std::vector<double> areas;
        STLVector<Eigen::Vector3d> point_cloud;
        double sum_area = 0.0;
        for(const Triangle& triangle : triangles)
        {
            double area = CalculateTriangleArea(triangle);
            areas.push_back(area);
            sum_area += area;
        }
        for(size_t i = 0; i < triangles.size(); i++)
        {
            if(visible_only)
            {
                Eigen::Vector3d dir = position - triangles[i].center;
                if(dir.dot(triangles[i].normal) < 0.0)
                    continue;
            }
            int num_points_i = std::ceil(areas[i] * num_points / sum_area);
            STLVector<Eigen::Vector3d> point_cloud_i = SampleTriangle(triangles[i], num_points_i);
            point_cloud.insert(point_cloud.end(),point_cloud_i.begin(),point_cloud_i.end());
        }
        return point_cloud;
    }

    static bool isPointInTriangle(const Eigen::Vector3d& point, const Triangle& triangle)
    {
        static const double eps = 1e-9;

        // if the point is on the plane
        Eigen::Vector3d AB = triangle.v2 - triangle.v1;
        Eigen::Vector3d AC = triangle.v3 - triangle.v1;
        Eigen::Vector3d AP = point - triangle.v1;
        Eigen::Vector3d normal = AB.cross(AC);
        if(fabs(normal.dot(AP)) > eps)
            return false;

        // if the point is in the triangle
        Eigen::Vector3d PA = - AP;
        Eigen::Vector3d PB = triangle.v2 - point;
        Eigen::Vector3d PC = triangle.v3 - point;

        Eigen::Vector3d n1 = PA.cross(PB);
        Eigen::Vector3d n2 = PB.cross(PC);
        Eigen::Vector3d n3 = PC.cross(PA);
        double dot_n1_n2 = n1.dot(n2);
        double dot_n2_n3 = n2.dot(n3);
        double dot_n3_n1 = n3.dot(n1);

        if(dot_n1_n2 * dot_n2_n3 > 0 && dot_n2_n3 * dot_n3_n1 > 0)
            return true;

        return false;
    }
    static STLVector<Eigen::Vector3d> SampleObjFile(const char* obj_file_name, int num_points)
    {
        STLVector<Eigen::Vector3d> vertices;
        ObjFileReader::ReadVerticesFromObj(obj_file_name, vertices);
        STLVector<Triangle> triangles;
        int n_triangles = vertices.size() / 3;
        for(int i=0; i<n_triangles; i++)
        {
            triangles.emplace_back(vertices[i*3],vertices[i*3+1],vertices[i*3+2]);
        }

        return SampleMultipleTriangles(triangles, num_points);
    }
};

}  // namespace generator
#endif