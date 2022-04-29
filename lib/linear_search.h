#ifndef _LINEAR_SEARCH_
#define _LINEAR_SEARCH_

#include "matrix.h"
#include <cmath>
#include <iostream>

double zoom(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const Matrix &initial, const Matrix & direction, 
            double a, double b, const double rho, const double sigma){
    double alpha, fa, ga;
    double f0 = f(initial), g0 = value(direction.T()*grad(initial));
    int step = 0;
    while(b-a>1e-8 && ++step<100){
        alpha = (a+b)/2;
        fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = value(direction.T()*grad(initial+alpha*direction));
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const Matrix &initial, const Matrix & direction, 
                    const double rho=0.4, const double sigma=0.7){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(initial), g0 = value(direction.T()*grad(initial));
    int step = 0;
    while(++step<100){
        double fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            return zoom(f, grad, initial, direction, a, alpha, rho, sigma);
        double ga = value(direction.T()*grad(initial+alpha*direction));
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom(f, grad, initial, direction, a, alpha, rho, sigma);
        a = alpha;
        alpha = (a+b)/2;
    }
    return alpha;
}

double armijo_goldstein(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), 
                        const Matrix &initial, const Matrix &direction, const double rho=0.4, const double t=100){
    double a = 0, b, alpha = 1;
    double f0 = f(initial), g0=value(direction.T()*grad(initial));
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

double simple_search(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const Matrix &initial, const Matrix & direction, 
                    const double rho, const double sigma){
    int step = 0;
    double alpha = 1, f0 = f(initial), g0 = value(direction.T()*grad(initial));
    while(++step < 50){
        double fa = f(initial + alpha * direction);
        if(fa <= f0 + sigma * alpha * g0)
            return alpha;
        alpha = alpha * rho;
    }
    return alpha;
}

double zoom_gradfree(double (*f)(const Matrix&), const Matrix &initial, const Matrix & direction, 
            double a, double b, const double rho, const double sigma){
    double alpha, fa, ga;
    double f0 = f(initial), g0 = value(direction.T()*gradient(f, initial));
    int step = 0;
    while(b-a>1e-8 && ++step<100){
        alpha = (a+b)/2;
        fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = value(direction.T()*gradient(f, initial+alpha*direction));
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell_gradfree(double (*f)(const Matrix&), const Matrix &initial, const Matrix & direction, 
                    const double rho=0.4, const double sigma=0.7){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(initial), g0 = value(direction.T()*gradient(f, initial));
    int step = 0;
    while(++step<100){
        double fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            return zoom_gradfree(f, initial, direction, a, alpha, rho, sigma);
        double ga = value(direction.T()*gradient(f, initial+alpha*direction));
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom_gradfree(f, initial, direction, a, alpha, rho, sigma);
        a = alpha;
        alpha = (a+b)/2;
    }
    return alpha;
}

double armijo_goldstein_gradfree(double (*f)(const Matrix&), const Matrix &initial, 
                        const Matrix &direction, const double rho=0.4, const double t=100){
    double a = 0, b, alpha = 1;
    double f0 = f(initial), g0=value(direction.T()*gradient(f,initial));
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