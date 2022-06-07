#ifndef _QUADRIC_PROGRAMMING_H_
#define _QUADRIC_PROGRAMMING_H_

#include "matrix.h"
#include <iostream>

/******************************************************************
 * 求解等式约束的凸二次规划问题，使用拉格朗日乘子法。问题表述：
 * min.  1/2 x^THx + c^Tx
 * s.t.  Ax=b
 * 调用方法：sol=quaprog_equ(H,c,A,b)
 * ****************************************************************/
Matrix quaprog_equ(const Matrix &H, const Matrix &c, const Matrix &A, const Matrix &b){
    Matrix IH = H.inv();
    Matrix AHA = A*IH*A.T();
    Matrix IAHA = AHA.inv();
    Matrix AIH = A*IH;
    Matrix G = IH-AIH.T()*IAHA*AIH;
    Matrix B = IAHA*AIH;
    return B.T()*b-G*c;
}

// 此函数用于求解有效集方法中的等式约束二次规划子问题，不建议单独调用
ColVector quaprog_subproblem(const Matrix &H, const ColVector &c, const Matrix &Ae, const ColVector &be, ColVector &lambda){
    Matrix ginvH = pinv(H);
    if(Ae.n){
        ColVector rb = Ae*ginvH*c+be;
        lambda = pinv(Ae*ginvH*Ae.T())*rb;
        return ginvH*(Ae.T()*lambda-c);
    } else {
        lambda = ColVector(be.n);
        return -ginvH*c;
    }
}

// 此函数用于统计p[l,...,r]中有多少个值为真
int checkednum(const bool *p, int l, int r){
    int ans = 0;
    for(int i = l; i <= r; i++)
        ans += p[i];
    return ans;
}

/******************************************************************
 * 求解等式约束的凸二次规划问题，使用有效集方法。问题表述：
 * min.  1/2 x^THx + c^Tx
 * s.t.  Ae x = be,  Ai x >= bi
 * 另有参数x0为指定的初始可行点
 * 调用方法：sol=quaprog(H,c,Ae,be,Ai,bi,x0,[MAXN],[eps])
 * ****************************************************************/
Matrix quaprog(const Matrix &H, const ColVector &c, const Matrix &Ae, const ColVector &be, 
               const Matrix &Ai, const ColVector &bi, const ColVector &x0, const int MAXN=1000, const double eps=1e-10){
    ColVector x = x0;
    int n = x.n, ne = be.n, ni = bi.n, step = 0;
    bool *idx = new bool[ni];
    for(int i = 0; i < ni; i++)
        idx[i] = ( fabs(Ai.getRow(i)*x-bi[i])<=eps );
    for(; step <= MAXN; step++){
        // 求解当前工作集对应的等式约束二次规划子问题
        Matrix Aee;
        if(ne) Aee = Ae;
        for(int j = 0; j < ni; j++)
            if(idx[j]) Aee = mergeRow(Aee,Ai.getRow(j));
        ColVector gk = H*x+c, lamk, dk;
        int m = Aee.n;
        dk = quaprog_subproblem(H,gk,Aee,zeroCol(m),lamk);
#ifdef DEBUG
        std::cerr << "step = " << step << std::endl;
        std::cerr << "x = " << x.T() << std::endl;
        std::cerr << "dk = " << dk.T() << std::endl;
        std::cerr << "lamk = " << lamk.T() << std::endl << std::endl;
#endif
        if(dk.vecnorm(2) <= eps){
            double y = 0.0; int jk;
            if(lamk.n > ne)
                for(int j = ne; j < lamk.n; j++)
                    if(lamk[j]<y) y=lamk[j], jk=j;  //寻找最小的乘子
            if(y >= 0) break;
            for(int i = 0; i < ni; i++)
                if(idx[i] && ne+checkednum(idx,0,i)==jk+1){  //将最小乘子对应的约束从工作集中移去
                    idx[i] = false;
                    break;
                }
        } else {
            int idxmin = -1;
            double alpha = 1.0;
            for(int i = 0; i < ni; i++){
                RowVector tmp = Ai.getRow(i);
                if( !idx[i] && tmp*dk<0 ){
                    double maxstep_i = (bi[i]-tmp*x)/(tmp*dk);
                    if(maxstep_i < alpha) alpha = maxstep_i, idxmin = i;  //寻找搜索方向上约束最紧的条件
                }
            }
            x = x + alpha*dk;
            if(idxmin != -1) idx[idxmin] = true;  //到达约束最紧的条件的边界，将其加入工作集
        }
    }
#ifndef SILENCE
    std::cerr << "---------- Quadric Programming Active Set Method ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << x.T() << std::endl;
    std::cerr << "OPtimal Value: " << 0.5*(x.T()*H*x) + c.T()*x << std::endl << std::endl;
#endif
    return x;
}

#endif