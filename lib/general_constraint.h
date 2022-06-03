#ifndef _GENERAL_CONSTRAINT_H_
#define _GENERAL_CONSTRAINT_H_

#include "matrix.h"
#include "quassi_newton.h"
#include <iostream>
#include <cmath>

namespace lagrange{
    static double (*f)(const ColVector&);
    static ColVector (*gradf)(const ColVector&);
    static ColVector (*h)(const ColVector&);
    static ColVector (*g)(const ColVector&);
    static Matrix (*Jh)(const ColVector&);
    static Matrix (*Jg)(const ColVector&);
    ColVector mu, lambda;
    double sigma;
    
    // 增广拉格朗日函数
    double psi(const ColVector &x){
        ColVector hval = h(x), gval = g(x);
        int l = hval.size(), m = gval.size();
        double res = f(x) - hval.T()*mu + 0.5*sigma*hval.sqrsum();
        double s2 = max(zeros(m,1),lambda-sigma*gval).sqrsum() - lambda.sqrsum();
        //std::cerr << "x=" << x.T() << "  psi=" << res + s2/(2.0*sigma) << std::endl;
        return res + s2/(2.0*sigma);
    }

    // 由于实现BFGS时技术落后，传参类型是Matrix而不是ColVector，故这里需要一个转化函数
    double psi(const Matrix &x){
        return psi((ColVector)x);
    }
    
    // 增广拉格朗日函数的梯度
    ColVector gradpsi(const ColVector &x){
        ColVector res=gradf(x);
        ColVector hval=h(x), gval=g(x);
        Matrix Jhval=Jh(x).T(), Jgval=Jg(x).T();
        int l=hval.size(), m=gval.size();
        for(int i = 0; i < l; i++)
            res = res + (sigma*hval[i]-mu[i])*Jhval.getCol(i);
        for(int i = 0; i < m; i++)
            res = res + (sigma*gval[i]-lambda[i])*Jgval.getCol(i);
        //std::cerr << "x=" << x.T() << "  grad=" << res.T() << std::endl;
        return res;
    }

    // 由于实现BFGS时技术落后，传参类型是Matrix而不是ColVector，故这里需要一个转化函数
    Matrix gradpsi(const Matrix &x){
        return gradpsi((ColVector)x);
    }
}

/*********************************************************************************
 * 求解一般约束优化问题的PHR算法通用子程序. 问题表述如下：
 *   min. f(x)
 *   s.t. h(x)=0, g(x)>=0
 * 其中h,g是向量函数，大于号是针对每个分量而言的
 * 
 * 使用方法： sol=general_constraint_optimal_PHR(f,h,g,gradf,Jh,Jg,x0,[eps],[theta],[eta],[sigma],[MAXN])
 * 必选参数：f,g,h含义见问题表述，gradf是f的梯度函数，Jh,Jg是h,g的Jacobi矩阵，x0是初始可行位置
 * 可选参数：eps是精度，theta和eta是PHR算法的实参数，sigma是罚因子，MAXN是最大迭代次数
 * *******************************************************************************/
ColVector general_constraint_optimal_PHR(double (*f)(const ColVector &x), ColVector (*h)(const ColVector &x), ColVector (*g)(const ColVector &x),
    ColVector (*gradf)(const ColVector &x), Matrix (*Jh)(const ColVector &x), Matrix (*Jg)(const ColVector &x), const ColVector &x0,
    const double eps=1e-6, const double theta=0.8, const double eta=2.0, const double sigma0=2.0, const int MAXN=5000){
    
    // 初始化函数指针
    lagrange::f = f;
    lagrange::h = h;
    lagrange::g = g;
    lagrange::gradf = gradf;
    lagrange::Jh = Jh;
    lagrange::Jg = Jg;
    lagrange::sigma = sigma0;
    using lagrange::mu;
    using lagrange::lambda;
    using lagrange::sigma;

    // 初始化迭代变量和乘子向量
    ColVector xold=x0, x;
    int n = x0.size();
    int l = h(x0).size();
    int m = g(x0).size();
    mu = 0.1*ones(l,1);
    lambda = 0.1*ones(m,1);
    double betak=10, betaold=10; // 这两个值用于检验终止条件
    int step=0;

    // 开始PHR算法的迭代过程
    for( ; betak>eps && step<MAXN; step++){
        x = bfgs(lagrange::psi, lagrange::gradpsi, xold, eps, 0.3, 0.7, 1000);
        ColVector hval=h(x), gval=g(x);
        betak = sqrt( hval.sqrsum() + min(gval,(1.0/sigma)*lambda).sqrsum() );
        if(betak>eps){
            mu  = mu-sigma*hval;
            lambda = max(zeroCol(m),lambda-sigma*gval);
            if(step>=2 && betak>theta*betaold) sigma *= eta;
        }
        betaold = betak;
        xold = x;
    }

#ifndef SILENCE
    std::cerr << "---------- General Constraint Optimal PHR Method ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << x.T() << std::endl;
    std::cerr << "OPtimal Value: " << f(x) << std::endl << std::endl;
#endif
    return x;
}

#endif