#pragma once
#ifndef UTILS_H_
#define UTILS_H_

#ifdef _WIN32
#include <io.h>
#elif __linux__
#include <unistd.h>
//#include <sys/types.h>
//#include <sys/stat.h>
#endif

// see https://www.cnblogs.com/renyu310/p/6485066.html for instructions of directory operation

class Utils
{
public:

    static void check_and_create_dir(const std::string& path)
    {
    #ifdef _WIN32
        if (_access(path.c_str(), 00) == -1)
        {
            std::string cmd = "mkdir " + unix_path_to_win(path);
            system(cmd.c_str());
        }
    #elif __linux__
        if (access(path.c_str(), 00) == -1)
        {
            std::string cmd = "mkdir -p " + path;
            system(cmd.c_str());
        }
    #endif
    }

private:

    #ifdef _WIN32
    static std::string unix_path_to_win(const std::string& unix_format)
    {
        std::string win_format(unix_format);
        for (int i = 0; i < win_format.length(); i++)
        {
            if (win_format[i] == '/')
            {
                win_format[i] = '\\';
            }
        }

        return win_format;
    }
    #endif

};

#endif  // UTILS_H_
