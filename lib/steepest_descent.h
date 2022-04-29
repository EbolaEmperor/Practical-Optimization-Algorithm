#ifndef _STEEPEST_DESCENT_H_
#define _STEEPEST_DESCENT_H_

#include <iostream>
#include "linear_search.h"
#include "matrix.h"

/*******************************************************************************************
 *
 * 这是一个最速下降法的通用最优化程序，带Wolfe搜索
 * 用法： x = steepest_descent(f,grad,hessian,[eps],[rho],[sigma])
 * 其中f是原函数，grad是梯度函数，hessian是Hessian矩阵函数，eps是收敛精度，rho和sigma是Wolfe准则的参数
 *
 *******************************************************************************************/
Matrix steepest_descent(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix (*hessian)(const Matrix&), Matrix current, 
    const double eps=1e-5, const double rho=0.1, const double sigma=0.4){
    int step = 0;
    while(grad(current).vecnorm(2)>eps)
    {
        step++;
        Matrix searchDirection = -grad(current);
        searchDirection = (1.0/searchDirection.vecnorm(2)) * searchDirection;
        double lambda = wolfe_powell(f,grad,current,searchDirection,rho,sigma,0.01*eps);
        current = current + lambda * searchDirection;
        //std::cout << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << std::endl;
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

#endif