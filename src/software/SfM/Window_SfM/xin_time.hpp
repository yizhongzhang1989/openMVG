//
// Created by xin on 2020/11/6.
//

#ifndef OPENMVG_XIN_TIME_HPP
#define OPENMVG_XIN_TIME_HPP

#include <ctime>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

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

    std::string time_str(  )
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        double second = elapsed_seconds.count();

        std::ostringstream out_put;
        out_put << "Cost Time: ";
//        std::cout << "Cost Time: ";
        if( second > 60 )
        {
            int64_t minute = static_cast<int64_t>(second / 60.);
            second = second - minute * 60;
            if( minute > 60 )
            {
                int64_t hour = static_cast<int64_t>(minute / 60.);
                minute = minute - hour * 60;
                out_put << hour << " hours ";
            }
            out_put << minute << " minutes ";
        }
        out_put << second << " seconds " << std::endl;

        return out_put.str();
    }

    void write_time( std::ofstream& out_file )
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        double second = elapsed_seconds.count();

        std::ostringstream out_put;
        out_put << "Cost Time: ";
//        std::cout << "Cost Time: ";
        if( second > 60 )
        {
            int64_t minute = static_cast<int64_t>(second / 60.);
            second = second - minute * 60;
            if( minute > 60 )
            {
                int64_t hour = static_cast<int64_t>(minute / 60.);
                minute = minute - hour * 60;
                out_put << hour << " hours ";
            }
            out_put << minute << " minutes ";
        }
        out_put << second << " seconds " << std::endl;

        std::cout << out_put.str() << std::endl;
//        out_put << out_put.rdbuf() ;
        out_file << out_put.str() ;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
};

#endif //OPENMVG_XIN_TIME_HPP
