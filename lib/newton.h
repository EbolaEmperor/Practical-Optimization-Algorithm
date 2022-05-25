#ifndef _NEWTON_H_
#define _NEWTON_H_

#include <iostream>
#include "linear_search.h"
#include "matrix.h"
#include "derivation.h"

/*********************************************************************************
 *
 * 这是一个纯牛顿法的通用最优化程序
 * 用法： x = newton(f,grad,hessian,x0,[eps])
 * 其中f是原函数，grad是梯度函数，hessian是Hessian矩阵函数，x0是初始位置，eps是收敛精度
 * 注意：牛顿法只适用于Hessian矩阵正定的函数
 *
 **********************************************************************************/
Matrix newton(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix (*hessian)(const Matrix&), Matrix current, const double eps=1e-5){
    int step = 0;
    while(grad(current).vecnorm(2)>eps){
        step++;
        Matrix searchDirection = solveByLDL(hessian(current), -grad(current));
        current = current + searchDirection;
        //std::cout << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << std::endl;
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

/*****************************************************************
 *
 * 这是一个纯牛顿法的通用最优化程序，无需人工求导
 * 用法： x = newton_gradfree(f,x0,[eps])
 * 其中f是原函数，x0是初始位置，eps是收敛精度
 * 注意：牛顿法只适用于Hessian矩阵正定的函数
 *
 ****************************************************************/
Matrix newton_gradfree(double (*f)(const Matrix&), Matrix current, const double eps=1e-5){
    int step = 0;
    while(1){
        step++;
        Matrix grad = gradient(f, current);
        if(grad.vecnorm(2)<=eps) break;
        Matrix searchDirection = solveByLDL(hesse(f,current), -grad);
        current = current + searchDirection;
#ifdef debug
        std::cerr << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << std::endl;
#endif
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

/*******************************************************************************************
 *
 * 这是一个带Wolfe不精确搜索的牛顿法通用最优化程序
 * 用法： x = newton_wolfe(f,grad,hessian,x0,[eps],[rho],[sigma])
 * 其中f是原函数，grad是梯度函数，hessian是Hessian矩阵函数，x0是初始位置，eps是收敛精度，rho和sigma是Wolfe准则的参数
 * 注意：牛顿法只适用于Hessian矩阵正定的函数
 *
 ******************************************************************************************/
Matrix newton_wolfe(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix (*hessian)(const Matrix&), Matrix current, 
    const double eps=1e-5, const double rho=0.1, const double sigma=0.4){
    int step = 0;
    while(grad(current).vecnorm(2)>eps){
        step++;
        Matrix searchDirection = solveByLDL(hessian(current), -grad(current));
        searchDirection = 1.0/searchDirection.vecnorm(2) * searchDirection;
        double lambda = wolfe_powell(f,grad,current,searchDirection,rho,sigma,eps);
        current = current + lambda * searchDirection;
        //std::cout << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << std::endl;
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

// Gill-Murray修正牛顿法的主程序
Matrix newton_gillmurray(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix (*hessian)(const Matrix&), Matrix current, 
    const double eps=1e-5, const double rho=0.1, const double sigma=0.4){
    int step = 0;
    while(grad(current).vecnorm(2)>eps){
        step++;
        Matrix searchDirection = solveByLDL_GM(hessian(current), -grad(current));
        searchDirection = 1.0/searchDirection.vecnorm(2) * searchDirection;
        double lambda = wolfe_powell(f,grad,current,searchDirection,rho,sigma,eps);
        current = current + searchDirection;
        std::cout << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << std::endl;
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    std::cout << "Total steps: " << step << std::endl;
    return current;
}

#endif