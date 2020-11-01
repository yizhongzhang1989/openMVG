//
// Created by v-xinli1 on 11/1/2020.
//

#ifndef OPENMVG_RANDOM_CIR_HPP
#define OPENMVG_RANDOM_CIR_HPP


class cRandom
{
public:
    cRandom(int x,double y):seed(x),random(y){};
    cRandom():seed(0),random(0){};

    int seed;
    double random;
};

cRandom my_random(int z)
// 16807 way to create random numbers
// z is the seed number, num is the total random number to create
{
    //z(n+1)=(a*z(n)+b) mod m
    //describe m=a*q+r to avoid that the number is large than the computer can bear
    const int m=std::pow(2,31)-1;
    const int a=16807;
    const int q=127773;
    const int r=2836;

    int temp=a*(z%q)-r*(z/q);

    if(temp<0)
    {
        temp=m+temp;
    }
    //z is the seed number
    z=temp;
    double t = z*1.0/m;

    cRandom cr;
    cr.random=t;
    cr.seed=z;

    return cr;
}



#endif //OPENMVG_RANDOM_CIR_HPP
