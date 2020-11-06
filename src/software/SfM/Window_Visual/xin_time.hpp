//
// Created by xin on 2020/11/6.
//

#ifndef OPENMVG_XIN_TIME_HPP
#define OPENMVG_XIN_TIME_HPP

#include <ctime>
#include <cstdlib>
#include <chrono>
#include <iostream>

class XinTime
{
public:
    XinTime()
    {
        tic();
    }

    void tic()
    {
        start = std::chrono::system_clock::now();
    }

    double toc()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() * 1000;
    }

    void print_time()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        double second = elapsed_seconds.count();
        std::cout << "Cost Time: ";
        if( second > 60 )
        {
            int64_t minute = static_cast<int64_t>(second / 60.);
            second = second - minute * 60;
            if( minute > 60 )
            {
                int64_t hour = static_cast<int64_t>(minute / 60.);
                minute = minute - hour * 60;
                std::cout << hour << " hours ";
            }
            std::cout << minute << " minutes ";
        }
        std::cout << second << " seconds " << std::endl;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
};

#endif //OPENMVG_XIN_TIME_HPP
