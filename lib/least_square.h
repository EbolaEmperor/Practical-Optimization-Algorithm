#ifndef _LEAST_SQUARE_
#define _LEAST_SQUARE_

#include "matrix.h"
#include "linear_search.h"
#include "derivation.h"
#include <cassert>

static ColVector (*curF)(const ColVector&);
static Matrix (*curJ)(const ColVector&);
static double curErr;

static double fval(const ColVector &x){
    ColVector y = curF(x);
    return y.T()*y/2;
}

static ColVector gradval(const ColVector &x){
    return curJ(x).T()*curF(x);
}

static ColVector gradval_gradfree(const ColVector &x){
    return jacobi(curF,x,0.05*curErr).T()*curF(x);
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
ColVector nonlinlsq_LM(ColVector (*F)(const ColVector&), Matrix (*J)(const ColVector&), ColVector current, const double err=1e-5, 
    const double rho=0.3, const double sigma=0.6, const int MAXN=2000){
    curF = F;
    curJ = J;
    const int n = current.n;
    int step = -1;
    double mu = fval(current);
    while(++step < MAXN){
        ColVector grad = gradval(current);
        if(grad.vecnorm(2) < err) break;
        Matrix Jk = J(current);
        ColVector Fk = F(current);
        // 搜索方向的确定是L-M方法的核心
        ColVector direct = -solve(Jk.T()*Jk + mu*eye(n), Jk.T()*Fk);
        double dfval = fval(current + direct) - fval(current);
        double q = Fk.T()*Jk*direct + 0.5*(Jk*direct).sqrsum();
        double r = dfval/q;
        // 若二次函数拟合效果较好，则减小mu以接近二次函数
        // 若二次函数拟合效果一般，则增大mu以限制dk的模长
        if(r>0.75) mu *= 0.1;
        else if(r<0.25) mu *= 10;
        // 计算调整mu后的搜索方向
        direct = -solve(Jk.T()*Jk + mu*eye(n), Jk.T()*Fk);
#ifdef DEBUG
        std::cerr << "step: " << step << "    current=" << current.T() << "    direction=" << direct.T() << "    f=" << fval(current) << std::endl;
#endif
        // 进行线搜索确定步长
        double alpha = wolfe_powell(fval, gradval, current, direct, rho, sigma, 0.05*err);
        current = current + alpha*direct;
    }
#ifndef SILENCE
    std::cerr << "---------- Least Square Problem Levenberg-Marquardt Method ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "OPtimal Value: " << fval(current) << std::endl << std::endl;
#endif
    return current;
}

/********************************************************************************************************
 *
 * 这是一个Levenberg-Marquardt方法求解非线性最小二乘问题的最优化程序，无需人工求导，一维搜索采用wolfe准则
 * 用法：sol=nonlinlsq_LM_gradfree(F,x0,[err],[rho],[sigma],[MAXN])
 * 其中F是待求解向量函数，x0是初始迭代位置，err是允许误差，
 * rho,sigma是Wolfe条件中的参数，MAXN是允许迭代的最大次数
 * 返回值是该方法求得的最小值点
 *
 ********************************************************************************************************/
ColVector nonlinlsq_LM_gradfree(ColVector (*F)(const ColVector&), ColVector current, const double err=1e-5, 
    const double rho=0.3, const double sigma=0.6, const int MAXN=2000){
    curF = F;
    curErr = err;
    const int n = current.n;
    int step = -1;
    double mu = fval(current);
    while(++step < MAXN){
        ColVector grad = gradval_gradfree(current);
        if(grad.vecnorm(2) < err) break;
        Matrix Jk = jacobi(F,current,0.05*err);
        ColVector Fk = F(current);
        ColVector direct = -solve(Jk.T()*Jk + mu*eye(n), Jk.T()*Fk);
#ifdef DEBUG
        std::cerr << "step: " << step << "    current=" << current.T() << "    direction=" << direct.T() << "    f=" << fval(current) << std::endl;
#endif
        double alpha = wolfe_powell(fval, gradval_gradfree, current, direct, rho, sigma, 0.05*err);
        double fval1 = fval(current);
        current = current + alpha*direct;
        double fval2 = fval(current);
        double q = Fk.T()*Jk*direct + 0.5*(Jk*direct).sqrsum();
        double r = (fval2-fval1)/q;
        if(r>0.75) mu *= 0.1;
        else if(r<0.25) mu *= 10;
    }
#ifndef SILENCE
    std::cerr << "---------- Least Square Problem Levenberg-Marquardt Method (gradfree) ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << current.T() << std::endl;
    std::cerr << "OPtimal Value: " << fval(current) << std::endl << std::endl;
#endif
    return current;
}

/********************************************************************************************************
 *
 * 这是一个 QR 分解法求解线性最小二乘问题的最优化程序
 * 用法：sol = linlsq(A, b)
 *
 ********************************************************************************************************/
ColVector linlsq(const Matrix &A, const ColVector &b){
    assert(A.n >= A.m);
    auto [Q, R] = A.getQR();
    int m = A.m;
    auto sol = solveUpperTriangular(R.getSubmatrix(0, m-1, 0, m-1), 
                                    (Q.T() * b).getSubmatrix(0, m-1, 0, 0));
#ifndef SILENCE
    std::cerr << "---------- Linear Least Square Problem ----------" << std::endl;
    std::cerr << "Optimal Point: " << sol.T() << std::endl;
    std::cerr << "OPtimal Value: " << norm(A * sol - b) << std::endl << std::endl;
#endif
    return sol;
}

#endif