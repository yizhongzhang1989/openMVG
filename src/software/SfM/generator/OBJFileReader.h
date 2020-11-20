#pragma once
#ifndef OBJ_FILE_READER_H_
#define OBJ_FILE_READER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#pragma warning(disable : 4996)

#include "types.h"

namespace generator
{

class ObjFileReader
{
public:
    static bool ReadVerticesFromObj(const char* obj_file_name, STLVector<Eigen::Vector3d>& vertices)
    {
        std::vector<double> v_list;
        if(readVertexListFromObj(obj_file_name,v_list))
        {
            if(v_list.size() % 3)
            {
                std::cerr<<"error reading: total numbers not divisible by 3."<<std::endl;
                return false;
            }
            int n_points = v_list.size() / 3;
            vertices.clear();
            for(int i = 0; i < n_points; i++)
            {
                vertices.emplace_back(v_list[i * 3], v_list[i * 3 + 1], v_list[i * 3 + 2]);
            }
        }
        else
            return false;

        return true;
    }

    static bool ReadVerticesOnlyFromObj(const char* obj_file_name, STLVector<Eigen::Vector3d>& vertices)
    {
        std::vector<double> v_list;
        if(readVertexOnlyFromObj(obj_file_name,v_list))
        {
            if(v_list.size() % 3)
            {
                std::cerr<<"error reading: total numbers not divisible by 3."<<std::endl;
                return false;
            }
            int n_points = v_list.size() / 3;
            vertices.clear();
            for(int i = 0; i < n_points; i++)
            {
                vertices.emplace_back(v_list[i * 3], v_list[i * 3 + 1], v_list[i * 3 + 2]);
            }
        }
        else
            return false;

        return true;
    }
private:
/**
* read obj file and create vertex list
* v_list format: f0_v0_x, f0_v0_y, f0_v0_z, f0_v1_x, f0_v1_y, f0_v1_z, f0_v2_x, f0_v2_y, f0_v2_z, f1_v0_x, f1_v0_y, f1_v0_z, ...
*/
    template<typename T>
    static inline int readVertexListFromObj(const char* obj_file_name, std::vector<T>& v_list)
    {
        std::ifstream obj(obj_file_name); //     load obj failed, do not touch old data
        if (!obj.is_open()) {
            std::cout << "error: readVertexListFromObj, cannot open " << obj_file_name << std::endl;
            return 0;
        }

        std::vector<T> v;
        std::vector<int> f;

        int face_format = 0; //     0:unkonwn;  001:v;  011:v/vt;  111:v/vt/vn;  101:v//vn;  1001:v//              1001 is not correct file format, but seen in some files
        std::string line;
        int line_idx = 0;
        while (std::getline(obj, line)) {
            line_idx++;
            std::stringstream ss;
            std::string cmd;
            ss << line;
            ss >> cmd;

            if (cmd == "v") {                 //     got a vertex, insert into v
                T xyz[3] = { 0, 0, 0 };
                int i = 0;
                while (i<3 && ss >> xyz[i])
                    i++;
                if (i < 3) {
                    std::cout << "read " << obj_file_name << " error:" << std::endl;
                    std::cout << "line " << line_idx << " : " << line << std::endl;
                    std::cout << "insert v anyway" << std::endl;
                }
                v.push_back(xyz[0]);
                v.push_back(xyz[1]);
                v.push_back(xyz[2]);
            }
            else if (cmd == "f") {            //     got a face, check its format, then insert into different parts
                std::string data;
                int fv_idx[3] = { 0, 0, 0 }, fvt_idx[3] = { 0, 0, 0 }, fvn_idx[3] = { 0, 0, 0 };
                int i = 0, count = 0;
                while (ss >> data) { //     we don't know how much vertices this face contain
                    if (face_format == 0) {                  //     face format unknown yet
                        face_format |= 0x01;       //     v
                        size_t found = data.find("//");
                        if (found != std::string::npos) {
                            if (found == data.length() - 2)
                                face_format |= 0x08; //     v//, an error case, but acceptable
                            else
                                face_format |= 0x04; //     v//vn
                        }
                        else {
                            found = data.find("/");
                            if (found != std::string::npos)
                                face_format |= 0x02; //     v/vt
                            found = data.find("/", found + 1);
                            if (found != std::string::npos)
                                face_format |= 0x04; //     v/vt/vn
                        }
                    }
                    if (i == 2) { //     more than 3 vertices, we need to clear the data
                        fv_idx[i] = fvt_idx[i] = fvn_idx[i] = 0;
                    }
                    int scan_check = 0;
                    if (face_format == 0x01 || face_format == 0x09) //     just v
                        scan_check = sscanf(data.c_str(), "%d", &fv_idx[i]);
                    else if (face_format == 0x03)                                 //     v/vt
                        scan_check = sscanf(data.c_str(), "%d/%d", &fv_idx[i], &fvt_idx[i]);
                    else if (face_format == 0x07)                                 //       v/vt/vn
                        scan_check = sscanf(data.c_str(), "%d/%d/%d", &fv_idx[i], &fvt_idx[i], &fvn_idx[i]);
                    else if (face_format == 0x05)                                 //     v//vn
                        scan_check = sscanf(data.c_str(), "%d//%d", &fv_idx[i], &fvn_idx[i]);
                    if (!scan_check) {
                        std::cout << "line " << line_idx << " : " << line << std::endl;
                        std::cout << "this vertex of this face is not read correctly due format or other error" << std::endl;
                        std::cout << "don't use this vertex" << std::endl;
                        face_format = 0;
                    }
                    else {
                        if (fv_idx[i] < 0)   //     negative index is relative index
                            fv_idx[i] = v.size() + fv_idx[i];
                        else
                            fv_idx[i] --; //     start from 1, so we have to minus 1 to make it start from 0

                        if (i == 2) {
                            if (face_format & 0x01) {
                                f.push_back(fv_idx[0]);
                                f.push_back(fv_idx[1]);
                                f.push_back(fv_idx[2]);
                            }
                            fv_idx[1] = fv_idx[2];
                            fvt_idx[1] = fvt_idx[2];
                            fvn_idx[1] = fvn_idx[2];
                            fv_idx[2] = 0;
                            fvt_idx[2] = 0;
                            fvn_idx[2] = 0;
                        }
                        else
                            i++;
                        count++;
                    }
                }
                if (count < 3) {
                    std::cout << "read " << obj_file_name << " error:" << std::endl;
                    std::cout << "line " << line_idx << " : " << line << std::endl;
                    std::cout << "too few vertices, don't insert f" << std::endl;
                }
            }
        }

        obj.close();

        //     create v_list
        v_list.clear();
        for (int i = 0; i < f.size(); i++) {
            int vid = f[i];
            v_list.push_back(v[vid * 3 + 0]);
            v_list.push_back(v[vid * 3 + 1]);
            v_list.push_back(v[vid * 3 + 2]);
        }

        return 1;
    }

    template<typename T>
    static inline bool readVertexOnlyFromObj(const char* obj_file_name, std::vector<T>& v_list)
    {
        std::ifstream obj(obj_file_name); // load obj failed, do not touch old data
        if (!obj.is_open())
        {
            std::cout << "error: readVertexListFromObj, cannot open " << obj_file_name << std::endl;
            return false;
        }

        std::vector<T> v;

        std::string line;
        int line_idx = 0;
        while (std::getline(obj, line))
        {
            line_idx++;
            std::stringstream ss;
            std::string cmd;
            ss << line;
            ss >> cmd;

            if (cmd == "v")
            {                 //     got a vertex, insert into v
                T xyz[3] = { 0, 0, 0 };
                int i = 0;
                while (i<3 && ss >> xyz[i])
                    i++;
                if (i < 3)
                {
                    std::cout << "read " << obj_file_name << " error:" << std::endl;
                    std::cout << "line " << line_idx << " : " << line << std::endl;
                    std::cout << "insert v anyway" << std::endl;
                }
                v.push_back(xyz[0]);
                v.push_back(xyz[1]);
                v.push_back(xyz[2]);
            }
        }

        obj.close();

        //  create v_list
        v_list = v;

        return true;
    }
};

}

#endif  // OBJ_FILE_READER_H_