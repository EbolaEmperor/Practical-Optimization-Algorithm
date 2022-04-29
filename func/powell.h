#ifndef _POWELL_H_
#define _POWELL_H_

#include "../lib/matrix.h"

namespace powell{

double f(const Matrix &x){
    return pow(x[0][0] + 10*x[1][0], 2)
           + 5 * pow(x[2][0] - x[3][0], 2)
           + pow(x[1][0] - 2*x[2][0], 4)
           + 10 * pow(x[0][0] - x[3][0], 4);
}

Matrix grad(const Matrix &x){
    Matrix g(4,1);
    g[0][0] = 2*(x[0][0]+10*x[1][0]) + 40*pow(x[0][0]-x[3][0],3);
    g[1][0] = 20*(x[0][0]+10*x[1][0]) + 4*pow(x[1][0]-2*x[2][0],3);
    g[2][0] = 10*(x[2][0]-x[3][0]) - 8*pow(x[1][0]-2*x[2][0],3);
    g[3][0] = -10*(x[2][0]-x[3][0]) - 40*pow(x[0][0]-x[3][0],3);
    return g;
}

Matrix hessian(const Matrix &x){
    Matrix h(4,4);
    h[0][0] = 2 + 120*pow(x[0][0]-x[3][0],2);
    h[0][1] = 20;
    h[0][3] = -120*pow(x[0][0]-x[3][0],2);
    h[1][0] = 20;
    h[1][1] = 200 + 12*pow(x[1][0]-2*x[2][0],2);
    h[1][2] = -24*pow(x[1][0]-2*x[2][0],2);
    h[2][1] = -24*pow(x[1][0]-2*x[2][0],2);
    h[2][2] = 10 + 48*pow(x[1][0]-2*x[2][0],2);
    h[2][3] = -10;
    h[3][0] = -120*pow(x[0][0]-x[3][0],2);
    h[3][2] = -10;
    h[3][3] = 10 + 120*pow(x[0][0]-x[3][0],2);
    return h;
}

Matrix initial(){
    static const double p[4] = {3,-1,0,1};
    return Matrix(p,4);
}

}

#endif