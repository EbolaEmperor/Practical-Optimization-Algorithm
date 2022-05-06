#ifndef _LEAST_SQUARE_
#define _LEAST_SQUARE_

#include "matrix.h"
#include "linear_search.h"

static Matrix (*curF)(const Matrix&);
static Matrix (*curJ)(const Matrix&);

static double fval(const Matrix &x){
    Matrix y = curF(x);
    return value(y.T()*y)/2;
}

static Matrix gradval(const Matrix &x){
    return curJ(x).T()*curF(x);
}

/********************************************************************************************************
 *
 * 这是一个Levenberg-Marquardt方法求解非线性最小二乘问题的最优化程序，一维搜索采用wolfe准则
 * 用法：sol=nonlinlsq_LM(F,J,x0,[err],[rho],[sigma],[MAXN])
 * 其中F是待求解向量函数，J是对应的Jacobi矩阵，x0是初始迭代位置，err是允许误差，
 * rho,sigma是Wolfe条件中的参数，MAXN是允许迭代的最大次数
 * 返回值是该方法求得的最小值点
 *
 ********************************************************************************************************/
Matrix nonlinlsq_LM(Matrix (*F)(const Matrix&), Matrix (*J)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.3, const double sigma=0.6, const int MAXN=2000){
    curF = F;
    curJ = J;
    const int n = current.n;
    int step = -1;
    double mu = fval(current);
    while(++step < MAXN){
        Matrix grad = gradval(current);
        if(grad.vecnorm(2) < err) break;
        Matrix Jk = J(current), Fk = F(current);
        Matrix direct = -solve(Jk.T()*Jk + mu*eye(n), Jk.T()*Fk);
        double alpha = wolfe_powell(fval, gradval, current, direct, rho, sigma, 0.05*err);
        double fval1 = fval(current);
        current = current + alpha*direct;
        double fval2 = fval(current);
        double q = value((Jk.T()*Fk).T()*direct) + 0.5*value((direct.T()*Jk.T()) * (Jk*direct));
        double r = (fval2-fval1)/q;
        if(r>0.75) mu *= 0.1;
        else if(r<0.25) mu *= 10;
    }
    std::cout << "Total Steps: " << step << std::endl;
    if(step==MAXN) std::cout << "Early Stop! The result might be unprecision." << std::endl;
    return current;
}

#endif