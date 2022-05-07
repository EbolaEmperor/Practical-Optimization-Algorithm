#ifndef _TRUST_REGION_
#define _TRUST_REGION_

#include "matrix.h"
#include "derivation.h"
#include <cmath>
#include <iostream>

Matrix dogleg(double (*f)(const Matrix&), const Matrix &g, const Matrix &B, const double &delta){
    Matrix pB = -solveByLDL(B,g);
    Matrix pU = -value(g.T()*g)/value(g.T()*B*g)*g;
    Matrix pBU = pB - pU;
    double tau;
    if(pB.vecnorm(2)<=delta)
        tau = 2;
    else if(pU.vecnorm(2)>delta)
        tau = delta/pU.vecnorm(2);
    else
        tau = sqrt(delta*delta-value(pU.T()*pU)) / value(pBU.T()*pBU) + 1;
    return (tau<=1) ? tau*pU : pU+(tau-1)*pBU;
}

Matrix trust_region_gradfree(double (*f)(const Matrix&), Matrix current, 
    const double deltamax=10, double delta=0.1, const double ita=0, const double err=1e-5){
    int step = 0;
    while(1){
        step++;
        Matrix g = gradient(f, current, 0.1*err);
        if(g.vecnorm(2)<err) break;
        Matrix B = hesse(f, current, 0.1*err);
        Matrix s = dogleg(f, g, B, delta);
        double rho = (f(current+s)-f(current))/(value(g.T()*s)+0.5*value(s.T()*B*s));
        if(rho < 0.25)
            delta *= 0.25;
        else if(rho > 0.75 && fabs(s.vecnorm(2)-delta)<1e-8)
            delta = std::min(2*delta, deltamax);
        if(rho > ita) current = current + s;
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

Matrix trust_region(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix (*hessian)(const Matrix&), Matrix current, 
    const double deltamax=10, double delta=0.1, const double ita=0, const double err=1e-5){
    int step = 0;
    while(1){
        step++;
        Matrix g = grad(current);
        if(g.vecnorm(2)<err) break;
        Matrix B = hessian(current);
        Matrix s = dogleg(f, g, B, delta);
        double rho = (f(current+s)-f(current))/(value(g.T()*s)+0.5*value(s.T()*B*s));
        if(rho < 0.25)
            delta *= 0.25;
        else if(rho > 0.75 && fabs(s.vecnorm(2)-delta)<1e-8)
            delta = std::min(2*delta, deltamax);
        if(rho > ita) current = current + s;
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

#endif