#ifndef _QUASSI_NEWTON_
#define _QUASSI_NEWTON_

#include "matrix.h"
#include "linear_search.h"
#include "derivation.h"
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
ColVector bfgs(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.3, const double sigma=0.6, const int MAXN=5000){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        if(step>MAXN) break;
        ColVector direction = -solveByLDL(B,grad(current));
        double norm = direction.vecnorm(2);
        direction = (1.0/norm)*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma,0.05*err);
        if(alpha < 0.1*err) alpha = norm;
        ColVector s = alpha*direction;
        ColVector y = grad(current+s) - grad(current);
        current = current + s;
        if(y.T()*s>0) //这是B保持正定的充要条件。事实上，因为一维搜索采用了Wolfe准则，这个条件总是成立，但其它准则下，这一条件的判断还是必要的
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << grad(current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal BFGS Method with Wolfe Condition ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector bfgs_goldstein(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.4, const int MAXN=5000){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        if(step>MAXN) break;
        ColVector direction = -solveByLDL(B,grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = armijo_goldstein(f,grad,current,direction,rho);
        ColVector s = alpha*direction;
        ColVector y = grad(current+s) - grad(current);
        current = current + s;
        if(y.T()*s>0)
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << grad(current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal BFGS Method with Goldstein Condition ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector bfgs_simple(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.55, const double sigma=0.4, const int MAXN=5000){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        if(step>MAXN) break;
        ColVector direction = -solveByLDL(B,grad(current));
        //direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = simple_search(f,grad,current,direction,rho,sigma);
        ColVector s = alpha*direction;
        ColVector y = grad(current+s) - grad(current);
        current = current + s;
        if(y.T()*s>0)
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << grad(current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal BFGS Method with Armijo Condition ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector dfp(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.55, const double sigma=0.4, const int MAXN=500000){
    const int n = current.n;
    Matrix H = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        if(step>MAXN) break;
        ColVector direction = -(H*grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma);
        ColVector s = alpha*direction;
        ColVector y = grad(current+s) - grad(current);
        current = current + s;
        H = H - (1.0/(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/(s.T()*s)) * (s*s.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << grad(current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal DFP Method with Wolfe Condition ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector broyden(double (*f)(const ColVector&), ColVector (*grad)(const ColVector&), ColVector current,
               const double phi=1, const double err=1e-5, const double rho=0.55, const double sigma=0.4, const int MAXN=50000){
    const int n = current.n;
    Matrix B = eye(n), H = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        ColVector direction = -solveByLDL(B,phi*grad(current)) - (1-phi)*(H*grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma);
        ColVector s = alpha*direction;
        ColVector y = grad(current+s) - grad(current);
        current = current + s;
        if(y.T()*s>0)
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
        H = H - (1.0/(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/(s.T()*s)) * (s*s.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << grad(current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal Broyden Method with Wolfe Condition ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector bfgs_gradfree(double (*f)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.3, const double sigma=0.6, const int MAXN=5000){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(gradient(f, current).vecnorm(2) > err){
        step++;
        ColVector direction = -solveByLDL(B,gradient(f, current));
        double norm = direction.vecnorm(2);
        direction = (1.0/norm)*direction;
        double alpha = wolfe_powell_gradfree(f,current,direction,rho,sigma);
        if(alpha < 0.1*err) alpha = norm;
        ColVector s = alpha*direction;
        ColVector y = gradient(f, current+s) - gradient(f, current);
        current = current + s;
        if(step > MAXN) break;
        if(y.T()*s>0)
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << gradient(f,current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal BFGS Method with Wolfe Condition (gradfree) ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector bfgs_goldstein_gradfree(double (*f)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.4, const int MAXN=5000){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        ColVector direction = -solveByLDL(B,gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = armijo_goldstein_gradfree(f,current,direction,rho);
        ColVector s = alpha*direction;
        ColVector y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        if(y.T()*s>0)
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << gradient(f,current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal BFGS Method with Goldstein Condition (gradfree) ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector dfp_gradfree(double (*f)(const ColVector&), ColVector current, const double err=1e-5, const double rho=0.55, const double sigma=0.4, const int MAXN=50000){
    const int n = current.n;
    Matrix H = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        ColVector direction = -(H*gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell_gradfree(f,current,direction,rho,sigma);
        ColVector s = alpha*direction;
        ColVector y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        H = H - (1.0/(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/(s.T()*s)) * (s*s.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << gradient(f,current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal DFP Method with Wolfe Condition (gradfree) ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
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
ColVector broyden_gradfree(double (*f)(const ColVector&), ColVector current, 
               const double phi=1, const double err=1e-5, const double rho=0.55, const double sigma=0.4, const int MAXN=50000){
    const int n = current.n;
    Matrix B = eye(n), H = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        ColVector direction = -solveByLDL(B,phi*gradient(f,current)) - (1-phi)*(H*gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell_gradfree(f,current,direction,rho,sigma);
        ColVector s = alpha*direction;
        ColVector y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        if(y.T()*s>0)
            B = B - (1.0/(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/(y.T()*s)) * (y*y.T());
        H = H - (1.0/(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/(s.T()*s)) * (s*s.T());
#ifdef DEBUG
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Step: " << step << std::endl;
        std::cout << "current = " << current.T() << std::endl;
        std::cout << "direction = " << direction.T() << std::endl;
        std::cout << "f = " << f(current)  << std::endl;
        std::cout << "||grad|| = " << gradient(f,current).vecnorm(2) << std::endl << std::endl;
#endif
    }
#ifndef SILENCE
    std::cerr << "---------- Non-Constraint Optimal Broyden Method with Wolfe Condition (gradfree) ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "Optimal Value: " << f(current) << std::endl << std::endl;
#endif
    return current;
}

#endif