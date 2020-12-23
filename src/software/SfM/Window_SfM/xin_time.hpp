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
    XinTime(bool start = true)
    {
        if (start)
            tic();
        elapsed_seconds_ = std::chrono::duration<double>(0);
    }

    void tic()
    {
        start_ = std::chrono::system_clock::now();
    }

    void toc()
    {
        end_ = std::chrono::system_clock::now();
        elapsed_seconds_ += (end_ - start_);
    }

    void pause()
    {
        end_ = std::chrono::system_clock::now();
        elapsed_seconds_ += (end_ - start_);
    }

    void restart()
    {
        start_ = std::chrono::system_clock::now();
    }

    void print_time()
    {
        double second = elapsed_seconds_.count();
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
        double second = elapsed_seconds_.count();

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
        double second = elapsed_seconds_.count();

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
    std::chrono::time_point<std::chrono::system_clock> start_, end_;
    std::chrono::duration<double> elapsed_seconds_;
};

#endif //OPENMVG_XIN_TIME_HPP
