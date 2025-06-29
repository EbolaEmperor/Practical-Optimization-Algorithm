#ifndef _LINEAR_SEARCH_
#define _LINEAR_SEARCH_

#include "matrix.h"
#include "derivation.h"
#include <cmath>
#include <iostream>

double zoom(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), const ColVector &initial, const ColVector & direction, 
            double a, double b, const double rho, const double sigma, const double eps=1e-8){
    double alpha, fa, ga;
    double f0 = f(initial), g0 = direction.T()*grad(initial);
    int step = 0;
    while(b-a>eps && ++step<100){
        alpha = (a+b)/2;
        fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = direction.T()*grad(initial+alpha*direction);
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), const ColVector &initial, const ColVector & direction, 
                    const double rho=0.4, const double sigma=0.7, const double eps=1e-8){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(initial), g0 = direction.T()*grad(initial);
    int step = 0;
    while(++step<100){
        double fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            return zoom(f, grad, initial, direction, a, alpha, rho, sigma, 0.1*eps);
        double ga = direction.T()*grad(initial+alpha*direction);
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom(f, grad, initial, direction, a, alpha, rho, sigma, 0.1*eps);
        a = alpha;
        alpha = (a+b)/2;
    }
    return alpha;
}

double armijo_goldstein(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), 
                        const ColVector &initial, const ColVector &direction, const double rho=0.4, const double t=100){
    double a = 0, b, alpha = 1;
    double f0 = f(initial), g0=direction.T()*grad(initial);
    bool b_is_inf = true;
    int maxstep = 200;
    while(1){
        double fa = f(initial + alpha*direction);
        //cout << a << " " << b << " " << alpha << " " << fa << " " << f0+rho*alpha*g0 << " " << f0+(1-rho)*alpha*g0 << endl;
        if(!--maxstep) return alpha;
        if(fa <= f0+rho*alpha*g0){
            if(fa >= f0+(1-rho)*alpha*g0) break;
            else{
                if(b_is_inf) alpha = t*alpha;
                else alpha = (a+b)/2;
            }
        } else {
            b = alpha;
            b_is_inf = false;
            alpha = (a+b)/2;
        }
    }
    return alpha;
}

double simple_search(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), const ColVector &initial, const ColVector &direction, 
                    const double rho, const double sigma){
    int step = 0;
    double alpha = 1, f0 = f(initial), g0 = direction.T()*grad(initial);
    while(++step < 20){
        double fa = f(initial + alpha * direction);
        if(fa <= f0 + sigma * alpha * g0)
            return alpha;
        alpha = alpha * rho;
    }
    return step==20 ? 1 : alpha;
}

double simple_search_gradfree(double (*f)(const ColVector&), const ColVector &initial, const ColVector &direction, 
                    const double rho, const double sigma, const double eps=1e-8){
    int step = 0;
    double alpha = 1, f0 = f(initial), g0 = direction_partial(f, initial, direction, eps);
    while(++step < 20){
        double fa = f(initial + alpha * direction);
        if(fa <= f0 + sigma * alpha * g0)
            return alpha;
        alpha = alpha * rho;
    }
    return step==20 ? 1 : alpha;
}

double zoom_gradfree(double (*f)(const ColVector&), const ColVector &initial, const ColVector &direction, 
            double a, double b, const double rho, const double sigma, const double eps){
    double alpha, fa, ga;
    double f0 = f(initial), g0 = direction_partial(f, initial, direction, eps);
    int step = 0;
    while(b-a>eps && ++step<100){
        alpha = (a+b)/2;
        fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = direction_partial(f, initial+alpha*direction, direction, eps);
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell_gradfree(double (*f)(const ColVector&), const ColVector &initial, const ColVector &direction, 
                    const double rho=0.4, const double sigma=0.7, const double eps = 1e-8){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(initial), g0 = direction_partial(f, initial, direction, eps);
    int step = 0;
    while(++step<100){
        double fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            return zoom_gradfree(f, initial, direction, a, alpha, rho, sigma, eps);
        double ga = direction_partial(f, initial+alpha*direction, direction, eps);
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom_gradfree(f, initial, direction, a, alpha, rho, sigma, eps);
        a = alpha;
        alpha = (a+b)/2;
    }
    return alpha;
}

double armijo_goldstein_gradfree(double (*f)(const ColVector&), const ColVector &initial, const ColVector &direction, 
                        const double rho=0.4, const double t=100){
    double a = 0, b, alpha = 1;
    double f0 = f(initial), g0=direction_partial(f, initial, direction);
    bool b_is_inf = true;
    int maxstep = 200;
    while(1){
        double fa = f(initial + alpha*direction);
        //cout << a << " " << b << " " << alpha << " " << fa << " " << f0+rho*alpha*g0 << " " << f0+(1-rho)*alpha*g0 << endl;
        if(!--maxstep) return alpha;
        if(fa <= f0+rho*alpha*g0){
            if(fa >= f0+(1-rho)*alpha*g0) break;
            else{
                if(b_is_inf) alpha = t*alpha;
                else alpha = (a+b)/2;
            }
        } else {
            b = alpha;
            b_is_inf = false;
            alpha = (a+b)/2;
        }
    }
    return alpha;
}

#endif