#ifndef _QUASSI_NEWTON_
#define _QUASSI_NEWTON_

#include "matrix.h"
#include "linear_search.h"
#include <cmath>
#include <iostream>

/********************************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，一维搜索采用wolfe准则
 * 用法：sol=bfgs(f,grad,x0,[err],[rho],[sigma])
 * 其中f是待求解函数，grad是待求解函数的梯度，x0是初始迭代位置，err是允许误差，rho,sigma是Wolfe条件中的参数
 * 返回值是BFGS方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix bfgs(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.3, const double sigma=0.6){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/********************************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，一维搜索采用goldstein准则
 * 用法：sol=bfgs_goldstein(f,grad,x0,[err],[rho])
 * 其中f是待求解函数，grad是待求解函数的梯度，x0是初始迭代位置，err是允许误差，rho是AG条件中的参数
 * 返回值是BFGS方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix bfgs_goldstein(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.4){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = armijo_goldstein(f,grad,current,direction,rho);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/********************************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，一维搜索采用简单准则
 * 用法：sol=bfgs_simple(f,grad,x0,[err],[rho],[sigma])
 * 其中f是待求解函数，grad是待求解函数的梯度，x0是初始迭代位置，err是允许误差，rho,sigma是简单准则中的参数
 * 返回值是BFGS方法求得的最小值点
 * 
 * 注：简单准则是指满足：f(xk + rho^m * dk) <= f(xk) + sigma * rho^m * gk^T * dk 的最小m
 * 找到m后令 x_k+1 = xk + rho^m *dk 进入下一次迭代
 *
 ********************************************************************************************************/
Matrix bfgs_simple(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = simple_search(f,grad,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/********************************************************************************************************
 *
 * 这是一个DFP拟牛顿法的通用最优化程序，一维搜索采用Wolfe准则
 * 用法：sol=dfp(f,grad,x0,[err],[rho],[sigma])
 * 其中f是待求解函数，grad是待求解函数的梯度，x0是初始迭代位置，err是允许误差，rho,sigma是简单准则中的参数
 * 返回值是DFP方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix dfp(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    const int n = current.n;
    Matrix H = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -(H*grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        H = H - (1.0/value(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/value(s.T()*s)) * (s*s.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/********************************************************************************************************
 *
 * 这是一个Broyden族校正拟牛顿法的通用最优化程序，一维搜索采用Wolfe准则
 * 用法：sol=broyden(f,grad,x0,[phi],[err],[rho],[sigma])
 * 其中f是待求解函数，grad是待求解函数的梯度，x0是初始迭代位置，err是允许误差，rho,sigma是Wolfe准则中的参数
 * phi是Broyden族校正中的重要参数，默认为1（即得BFGS校正），为0时就是DFP校正
 * 返回值是Broyden族校正方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix broyden(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), Matrix current, 
               const double phi=1, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    const int n = current.n;
    Matrix B = eye(n), H = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,phi*grad(current)) - (1-phi)*(H*grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
        H = H - (1.0/value(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/value(s.T()*y)) * (s*s.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/***********************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，可自动求导，一维搜索采用Wolfe准则
 * 用法：sol=bfgs(f,x0,[err],[rho],[sigma])
 * 其中f是待求解函数,x0是初始迭代位置，err是允许误差，rho,sigma是Wolfe条件中的参数
 * 返回值是BFGS方法求得的最小值点
 *
 ***********************************************************************************************/
Matrix bfgs_gradfree(double (*f)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.3, const double sigma=0.6){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(gradient(f, current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,gradient(f, current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell_gradfree(f,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = gradient(f, current+s) - gradient(f, current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/*****************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，可自动求导，一维搜索采用Goldstein准则
 * 用法：sol=bfgs(f,x0,[err],[rho])
 * 其中f是待求解函数，x0是初始迭代位置，err是允许误差，rho是AG条件中的参数
 * 返回值是BFGS方法求得的最小值点
 *
 ****************************************************************************************/
Matrix bfgs_goldstein_gradfree(double (*f)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.4){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = armijo_goldstein_gradfree(f,current,direction,rho);
        Matrix s = alpha*direction;
        Matrix y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/********************************************************************************************************
 *
 * 这是一个DFP拟牛顿法的通用最优化程序，可自动求导，一维搜索采用Wolfe准则
 * 用法：sol=dfp(f,x0,[err],[rho],[sigma])
 * 其中f是待求解函数，x0是初始迭代位置，err是允许误差，rho,sigma是简单准则中的参数
 * 返回值是DFP方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix dfp_gradfree(double (*f)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    const int n = current.n;
    Matrix H = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        Matrix direction = -(H*gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell_gradfree(f,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        H = H - (1.0/value(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/value(s.T()*s)) * (s*s.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

/********************************************************************************************************
 *
 * 这是一个Broyden族校正拟牛顿法的通用最优化程序，可自动求导，一维搜索采用Wolfe准则
 * 用法：sol=broyden(f,x0,[phi],[err],[rho],[sigma])
 * 其中f是待求解函数，x0是初始迭代位置，err是允许误差，rho,sigma是Wolfe准则中的参数
 * phi是Broyden族校正中的重要参数，默认为1（即得BFGS校正），为0时就是DFP校正
 * 返回值是Broyden族校正方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix broyden_gradfree(double (*f)(const Matrix&), Matrix current, 
               const double phi=1, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    const int n = current.n;
    Matrix B = eye(n), H = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,phi*gradient(f,current)) - (1-phi)*(H*gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell_gradfree(f,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << gradient(f,current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
        H = H - (1.0/value(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/value(s.T()*y)) * (s*s.T());
    }
    std::cout << "Total Steps: " << step << std::endl;
    return current;
}

#endif