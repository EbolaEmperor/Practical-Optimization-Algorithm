#ifndef _CONGUGATE_GRADIENT_
#define _CONGUGATE_GRADIENT_

#include "matrix.h"
#include "linear_search.h"
#include "derivation.h"
#include "preconditioner.h"
#include <iostream>

/***************************************************************
 *
 * 这是一个用朴素共轭梯度法求解方程Ax=b的通用最优化程序
 * 调用CG(A,b,x,err)即可，其中x为初始迭代位置
 *
 **************************************************************/
Matrix CG(const Matrix &A, const Matrix &b, Matrix x, const double err = 1e-6){
    Matrix r = A*x-b;
    Matrix p = -r;
    long long step = 0;
    while(r.vecnorm(2) >= err){
        step++;
        double alpha = - value(r.T()*p) / value(p.T()*A*p);
        x = x + alpha*p;
        double tmp = value(r.T()*r);
        r = r + alpha*A*p;
        double beta = value(r.T()*r) / tmp;
        p = -r + beta*p;
    }
    std::cout << "Steps: " << step << std::endl;
    return x;
}

/***************************************************************
 *
 * 这是一个用预优共轭梯度法求解方程Ax=b的通用最优化程序
 * 调用PCG(A,b,x,M,err)即可，其中x为初始迭代位置, M是预优因子
 *
 **************************************************************/
template <typename Mat>
ColVector PCG(const Mat &A, const ColVector &b, ColVector x, const Preconditioner &P, const double err = 1e-6){
    ColVector r = A * x - b, y = P.vmult(r), p = -y;
    long long step = 0;
    double dotry = dot(r, y);
    while(r.vecnorm(2) >= err){
        step++;
        auto Ap = A * p;
        double alpha = dotry / dot(p, Ap);
        x += alpha * p;
        r += alpha * Ap;
        y = P.vmult(r);
        double tmp = dotry;
        dotry = dot(r, y);
        double beta = dotry / tmp;
        (p *= beta) -= y;
    }
    std::cout << "PCG Steps: " << step << std::endl;
    return x;
}

/**************************************************************************
 *
 * 这是一个用FR非线性共轭梯度法的通用最优化程序
 * 调用CG_FR(f,grad,x,err)即可，其中f为待优化函数，grad为梯度，x为初始迭代位置
 *
 **************************************************************************/
ColVector CG_FR(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), ColVector x, const double err = 1e-6, const double rho=0.3, const double sigma=0.6){
    ColVector g = grad(x);
    ColVector d = -g;
    int step = 0;
    const int n = x.n;
    while(g.vecnorm(2) >= err){
        step++;
        double alpha = wolfe_powell(f, grad, x, d, rho, sigma, 0.05*err);
        x = x + alpha*d;
        if(step % (n+1) == 0){
            g = grad(x);
            d = -g;
        } else {
            double tmp = g.sqrsum();
            g = grad(x);
            double beta = g.sqrsum() / tmp;
            d = -g + beta * d;
        }
    }
    std::cout << "Steps: " << step << std::endl;
    return x;
}

/**************************************************************************
 *
 * 这是一个用FR非线性共轭梯度法的通用最优化程序，无需人工求导
 * 调用CG_FR(f,grad,x,err)即可，其中f为待优化函数，grad为梯度，x为初始迭代位置
 *
 **************************************************************************/
ColVector CG_FR_gradfree(double (*f)(const ColVector&), ColVector x, const double err = 1e-6, const double rho=0.3, const double sigma=0.6){
    ColVector g = gradient(f, x);
    ColVector d = -g;
    int step = 0;
    const int n = x.n;
    while(g.vecnorm(2) >= err){
        step++;
        double alpha = wolfe_powell_gradfree(f, x, d, rho, sigma, 0.05*err);
        x = x + alpha*d;
        if(step % (n+1) == 0){
            g = gradient(f, x);
            d = -g;
        } else {
            double tmp = g.sqrsum();
            g = gradient(f, x);
            double beta = g.sqrsum() / tmp;
            d = -g + beta * d;
        }
    }
    std::cout << "Steps: " << step << std::endl;
    return x;
}

#endif