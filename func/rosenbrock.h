#ifndef _ROSENBROCK_
#define _ROSENBROCK_

#include "../lib/matrix.h"

namespace rosenbrock{

int n;

double f(const ColVector &x){
    double res = 0;
    for(int i = 0; i < n-1; i++)
        res += 100*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (1-x[i])*(1-x[i]);
    return res;
}

ColVector grad(const ColVector &x){
    ColVector g(n);
    for(int i = 0; i < n-1; i++){
        g[i]   += 400*x[i]*(x[i]*x[i]-x[i+1]) + 2*(x[i]-1);
        g[i+1] += 200*(x[i+1]-x[i]*x[i]);
    }
    return g;
}

Matrix hessian(const ColVector &x){
    Matrix h(n,n);
    h[0][0] = 1200*x[0]*x[0] - 400*x[1] + 2;
    h[0][1] = -400*x[0];
    for(int i = 1; i < n-1; i++){
        h[i][i-1] = -400*x[i-1];
        h[i][i] = 1200*x[i]*x[i] - 400*x[i+1] + 202;
        h[i][i+1] = -400*x[i];
    }
    h[n-1][n-2] = -400*x[n-2];
    h[n-1][n-1] = 200;
    return h;
}

ColVector initial(){
    ColVector x(n);
    for(int i = 0; i < n; i += 2){
        x[i] = -1.2;
        if(i+1<n) x[i+1] = 1;
    }
    return x;
}

}

#endif