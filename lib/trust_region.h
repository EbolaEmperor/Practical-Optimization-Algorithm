#ifndef _TRUST_REGION_
#define _TRUST_REGION_

#include "matrix.h"
#include "derivation.h"
#include <cmath>
#include <iostream>

ColVector dogleg(double (*f)(const ColVector&), const ColVector &g, const Matrix &B, const double &delta){
    ColVector pB = -solveByLDL(B,g);
    ColVector pU = -(g.T()*g)/(g.T()*B*g) * g;
    ColVector pBU = pB - pU;
    double tau;
    if(pB.vecnorm(2)<=delta)
        tau = 2;
    else if(pU.vecnorm(2)>delta)
        tau = delta/pU.vecnorm(2);
    else
        tau = sqrt(delta*delta-pU.sqrsum()) / pBU.sqrsum() + 1;
    return (tau<=1) ? tau*pU : pU+(tau-1)*pBU;
}

ColVector trust_region_gradfree(double (*f)(const ColVector&), ColVector current, 
    const double deltamax=10, double delta=0.1, const double ita=0, const double err=1e-5, const int MAXN=10000){
    int step = 0;
    while(1){
        step++;
        ColVector g = gradient(f, current, 0.1*err);
        if(g.vecnorm(2)<err) break;
        Matrix B = hesse(f, current, 0.1*err);
        ColVector s = dogleg(f, g, B, delta);
        double rho = (f(current+s)-f(current))/(g.T()*s+0.5*(s.T()*B*s));
        if(rho < 0.25)
            delta *= 0.25;
        else if(rho > 0.75 && fabs(s.vecnorm(2)-delta)<1e-8)
            delta = std::min(2*delta, deltamax);
        if(rho > ita) current = current + s;
    }

#ifndef SILENCE
    std::cerr << "---------- Non-constraint Optimal Trust Region Method (gradfree) ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "OPtimal Value: " << f(current) << std::endl << std::endl;
#endif

    return current;
}

ColVector trust_region(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), Matrix (*hessian)(const ColVector&), ColVector current, 
    const double deltamax=10, double delta=0.1, const double ita=0, const double err=1e-5, const int MAXN=10000){
    int step = 0;
    while(1){
        step++;
        ColVector g = grad(current);
        if(g.vecnorm(2)<err) break;
        Matrix B = hessian(current);
        ColVector s = dogleg(f, g, B, delta);
        double rho = (f(current+s)-f(current))/(g.T()*s+0.5*(s.T()*B*s));
        if(rho < 0.25)
            delta *= 0.25;
        else if(rho > 0.75 && fabs(s.vecnorm(2)-delta)<1e-8)
            delta = std::min(2*delta, deltamax);
        if(rho > ita) current = current + s;
    }

#ifndef SILENCE
    std::cerr << "---------- Non-constraint Optimal Trust Region Method ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "OPtimal Value: " << f(current) << std::endl << std::endl;
#endif

    return current;
}

#endif