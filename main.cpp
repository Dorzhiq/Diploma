#include <iostream>
#include <math.h>

double fun(double x, double y, double z){
    return x+y+z;
}

double dif2(double x, double y, double z, double(*U)(double,double,double),
            double h1,double h2, double h3){
    return (U(x,y,z) + U(x+h1, y+h2, z+h3))/2;
}

int main(int argc, char *argv[])
{
    double W[10][10];
    std::cout<<dif2(1,1,1,fun,0.1,0.05,0.01)<<std::endl;
    std::cout<<dif2(2,2,2,fun,0.1,0.05,0.01)<<std::endl;
}

