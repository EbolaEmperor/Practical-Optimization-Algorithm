#ifndef _POWELL_H_
#define _POWELL_H_

#include "../lib/matrix.h"

namespace powell{

double f(const ColVector &x){
    return pow(x[0] + 10*x[1], 2)
           + 5 * pow(x[2] - x[3], 2)
           + pow(x[1] - 2*x[2], 4)
           + 10 * pow(x[0] - x[3], 4);
}

ColVector grad(const ColVector & x){
    ColVector g(4);
    g[0] = 2*(x[0]+10*x[1]) + 40*pow(x[0]-x[3],3);
    g[1] = 20*(x[0]+10*x[1]) + 4*pow(x[1]-2*x[2],3);
    g[2] = 10*(x[2]-x[3]) - 8*pow(x[1]-2*x[2],3);
    g[3] = -10*(x[2]-x[3]) - 40*pow(x[0]-x[3],3);
    return g;
}

Matrix hessian(const ColVector &x){
    Matrix h(4,4);
    h[0][0] = 2 + 120*pow(x[0]-x[3],2);
    h[0][1] = 20;
    h[0][3] = -120*pow(x[0]-x[3],2);
    h[1][0] = 20;
    h[1][1] = 200 + 12*pow(x[1]-2*x[2],2);
    h[1][2] = -24*pow(x[1]-2*x[2],2);
    h[2][1] = -24*pow(x[1]-2*x[2],2);
    h[2][2] = 10 + 48*pow(x[1]-2*x[2],2);
    h[2][3] = -10;
    h[3][0] = -120*pow(x[0]-x[3],2);
    h[3][2] = -10;
    h[3][3] = 10 + 120*pow(x[0]-x[3],2);
    return h;
}

ColVector initial(){
    static const double p[4] = {3,-1,0,1};
    return ColVector(4,p);
}

}

#endif