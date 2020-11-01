#pragma once
#ifndef COLOR_GENERATOR_H_
#define COLOR_GENERATOR_H_

#include <random>
#include <chrono>
#include "types.h"

namespace generator
{

template<typename ColorType>
class ColorGenerator
{
public:
    ColorGenerator()
            : e(std::chrono::system_clock::now().time_since_epoch().count())
    {
        n = std::uniform_int_distribution<unsigned char>(0,255);
    }
    ColorType Generate()
    {
        unsigned char r = n(e);
        unsigned char g = n(e);
        unsigned char b = n(e);
        return {r,g,b};
    }
private:
    std::default_random_engine e;
    std::uniform_int_distribution<unsigned char> n;
};

}  // namespace generator


#endif  // COLOR_GENERATOR_H_