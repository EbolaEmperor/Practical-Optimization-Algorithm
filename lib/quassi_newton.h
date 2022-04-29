#ifndef _QUASSI_NEWTON_
#define _QUASSI_NEWTON_

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

// 分割线---以下是gradfree版本
//***************************************************************************************************

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