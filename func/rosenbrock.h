#ifndef _ROSENBROCK_
#define _ROSENBROCK_

#include "../lib/matrix.h"

namespace rosenbrock{

int n;

double f(const Matrix &x){
    double res = 0;
    for(int i = 0; i < n-1; i++)
        res += 100*(x[i+1][0]-x[i][0]*x[i][0])*(x[i+1][0]-x[i][0]*x[i][0]) + (1-x[i][0])*(1-x[i][0]);
    return res;
}

Matrix grad(const Matrix &x){
    Matrix g(n,1);
    for(int i = 0; i < n-1; i++){
        g[i][0]   += 400*x[i][0]*(x[i][0]*x[i][0]-x[i+1][0]) + 2*(x[i][0]-1);
        g[i+1][0] += 200*(x[i+1][0]-x[i][0]*x[i][0]);
    }
    return g;
}

Matrix hessian(const Matrix &x){
    Matrix h(n,n);
    h[0][0] = 1200*x[0][0]*x[0][0] - 400*x[1][0] + 2;
    h[0][1] = -400*x[0][0];
    for(int i = 1; i < n-1; i++){
        h[i][i-1] = -400*x[i-1][0];
        h[i][i] = 1200*x[i][0]*x[i][0] - 400*x[i+1][0] + 202;
        h[i][i+1] = -400*x[i][0];
    }
    h[n-1][n-2] = -400*x[n-2][0];
    h[n-1][n-1] = 200;
    return h;
}

}

#endif