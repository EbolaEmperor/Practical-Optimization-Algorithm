#ifndef _STEEPEST_DESCENT_H_
#define _STEEPEST_DESCENT_H_

#include <iostream>
#include "linear_search.h"
#include "findzero.h"
#include "matrix.h"

/*******************************************************************************************
 *
 * 这是一个最速下降法的通用最优化程序，带Wolfe搜索
 * 用法： x = steepest_descent(f,grad,hessian,[eps],[rho],[sigma])
 * 其中f是原函数，grad是梯度函数，hessian是Hessian矩阵函数，eps是收敛精度，rho和sigma是Wolfe准则的参数
 *
 *******************************************************************************************/
Matrix steepest_descent(double (*f)(const ColVector&), 
                        ColVector (*grad)(const ColVector&), 
                        ColVector current,
                        const double eps=1e-5, 
                        const double rho=0.1, 
                        const double sigma=0.4){
    int step = 0;
    while(grad(current).vecnorm(2)>eps)
    {
        step++;
        ColVector searchDirection = -grad(current);
        searchDirection = (1.0/searchDirection.vecnorm(2)) * searchDirection;
        double lambda = wolfe_powell(f,grad,current,searchDirection,rho,sigma,0.01*eps);
        current = current + lambda * searchDirection;
        //std::cout << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << std::endl;
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

/*******************************************************************************************
 *
 * 这是一个最速下降法的通用最优化程序，一维搜索使用割线法，对方向导函数的零点进行精确搜索
 * 用法： x = steepest_descent_exact(f,grad,hessian,[eps])
 * 其中f是原函数，grad是梯度函数，hessian是Hessian矩阵函数，eps是收敛精度
 *
 *******************************************************************************************/
Matrix steepest_descent_exact(double (*f)(const ColVector&), 
                              ColVector (*grad)(const ColVector&), 
                              ColVector current, 
                              const double eps=1e-5){
    int step = 0;
    while(grad(current).vecnorm(2)>eps)
    {
        step++;
        ColVector searchDirection = -grad(current);
        searchDirection = (1.0/searchDirection.vecnorm(2)) * searchDirection;

        auto f_1dim = [&](double alpha) {
            // f1(a) = f(x + ad),  f1'(a) = d' * grad(x + ad)
            return dot(searchDirection, grad(current + alpha * searchDirection));
        };
        double alpha = findzero_secant(f_1dim, 0.0, 1.0, 0.1*eps);

        current = current + alpha * searchDirection;
#ifdef DEBUG
        std::cout << "step: " << step << "    x=" << current.T() << "    direction=" 
                  << searchDirection.T() << "     alpha=" << alpha << "    f(x)=" << f(current) << std::endl;
#endif
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

#endif