#ifndef _DERIVATION_H_
#define _DERIVATION_H_

#include <iostream>
#include <cmath>
#include "matrix.h"

double partial(double (*f)(const ColVector&), const ColVector &x, const int &j, const double err=1e-6){
    const int n = x.n;
    if(j<0 || j>=n){
        std::cerr << "partial():: out of range" << std::endl;
        return 0;
    }
    ColVector d(n);
    d[j] = 1;
    double alpha = 0.5;
    double g1 = alpha*(f(x+d)-f(x-d)), g2 = g1+10*err;
    int step = 0;
    while(fabs(g2-g1)>err && ++step<12){
        alpha *= 10;
        d = 0.1 * d;
        g2 = g1;
        g1 = alpha*(f(x+d)-f(x-d));
    }
    return g1;
}

double direction_partial(double (*f)(const ColVector&), const ColVector &x, ColVector d, const double err=1e-6){
    const int n = x.n;
    if(d.n != n){
        std::cerr << "direction_partial():: incompatible dimensions" << std::endl;
        return 0;
    }
    d = d / d.vecnorm(2) * (err * 10);
    double alpha = 0.5 / (err * 10);
    double g1 = alpha*(f(x+d)-f(x-d)), g2 = g1+10*err;
    int step = 0;
    while(fabs(g2-g1)>err && ++step<12){
        alpha *= 10;
        d = 0.1 * d;
        g2 = g1;
        g1 = alpha*(f(x+d)-f(x-d));
    }
    return g1;
}

ColVector gradient(double (*f)(const ColVector&), const ColVector &x, const double err=1e-6){
    const int n = x.n;
    ColVector g(n);
    for(int i = 0; i < n; i++)
        g[i] = partial(f, x, i, err);
    return g;
}

double partial2(double (*f)(const ColVector&), const ColVector &x, const int &i, const int &j, const double err=1e-6){
    const int n = x.n;
    if(j<0 || j>=n || i<0 || i>=n){
        std::cerr << "partial():: out of range" << std::endl;
        return 0;
    }
    ColVector d1(n), d2(n);
    d1[i] = 1;
    d2[j] = 1;
    double alpha = 0.5;
    double g1 = alpha*alpha*(f(x+d1+d2)+f(x-d1-d2)-f(x+d1-d2)-f(x-d1+d2)), g2 = g1+10*err;
    int step = 0;
    while(fabs(g2-g1)>err && ++step<5){
        alpha *= 10;
        d1 = 0.1 * d1;
        d2 = 0.1 * d2;
        g2 = g1;
        g1 = alpha*alpha*(f(x+d1+d2)+f(x-d1-d2)-f(x+d1-d2)-f(x-d1+d2));
    }
    return g1;
}

Matrix hesse(double (*f)(const ColVector&), const ColVector &x, const double err=1e-6){
    const int n = x.n;
    Matrix h(n,n);
    for(int i = 0; i < n; i++){
        h[i][i] = partial2(f, x, i, i, err);
        for(int j = 0; j < i; j++)
            h[j][i] = h[i][j] = partial2(f, x, i, j, err);
    }
    return h;
}

Matrix jacobi(ColVector (*F)(const ColVector&), const ColVector &x, const double err=1e-6){
    ColVector y0 = F(x);
    const int n = x.n;
    const int m = y0.n;
    Matrix h(m,n);
    for(int j = 0; j < n; j++){
        ColVector d(n);
        d[j] = 1;
        double alpha = 0.5;
        ColVector g1 = alpha*(F(x+d)-F(x-d));
        double g2 = g1.vecnorm(2)+10*err;
        int step = 0;
        while(fabs(g2-g1.vecnorm(2))>err && ++step<12){
            alpha *= 10;
            d = 0.1 * d;
            g2 = g1.vecnorm(2);
            g1 = alpha*(F(x+d)-F(x-d));
        }
        h.setSubmatrix(0,j,g1);
    }
    return h;
}

#endif