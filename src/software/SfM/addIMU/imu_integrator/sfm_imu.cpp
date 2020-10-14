// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sfm_imu.hpp"
#include <fstream>
#include <sstream>

namespace openMVG {
namespace sfm {

class csvReader
{
public:
    csvReader(const std::string& filename);
    ~csvReader();
    bool readline();
    float data[10];
private:
    std::ifstream _csvInput;

};  // csvReader

csvReader::csvReader(const std::string& filename)
{
    _csvInput.open(filename.c_str());
    if(!_csvInput)
        throw std::runtime_error("Can not open csv file "+filename+"\n");
    std::string header;
    getline(_csvInput,header);
}

csvReader::~csvReader()
{
    if(_csvInput)
        _csvInput.close();
}

bool csvReader::readline()
{
    if(!_csvInput)
        return false;
    if(_csvInput.eof())
        return false;

    std::string line;
    getline(_csvInput,line);

    std::istringstream _readStr(line);
    std::string part;


    for(int i=0;i<10;i++)
    {
        getline(_readStr,part,',');
        data[i] = atof(part.c_str());
    }

//    for(int i=0;i<7;i++)
//    {
//        getline(_readStr,part,',');
//        data[i] = atof(part.c_str());
//    }
//    for(int i=7;i<10;i++)
//    {
//        data[i] = 0;
//    }
//    for(int i=1;i<4;i++)
//    {
//        double temp = data[i];
//        data[i] = data[i+3];
//        data[i+3] = temp;
//    }

    bool flag = false;
    for(int i=1;i<10;i++)
    {
        if(data[i])
        {
            flag = true;
            break;
        }
    }

    return flag;
}

void SfM_IMU::ReadFromCSV(const std::string& filename)
{
    csvReader reader(filename);

    while(reader.readline())
    {
        timestamps.push_back(reader.data[0]);
        vec_acc.emplace_back(reader.data[1],reader.data[2],reader.data[3]);
        vec_gyro.emplace_back(reader.data[4],reader.data[5],reader.data[6]);
    }
}


}  // namespace sfm
}  // namespace openMVG