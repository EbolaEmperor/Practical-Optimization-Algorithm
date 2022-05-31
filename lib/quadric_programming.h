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

Matrix quaprog_subproblem(const Matrix &H, const Matrix &c, const Matrix &Ae, const Matrix &be, Matrix &lambda){
    Matrix ginvH = pinv(H);
    if(Ae.n){
        Matrix rb = Ae*ginvH*c+be;
        lambda = pinv(Ae*ginvH*Ae.T())*rb;
        return ginvH*(Ae.T()*lambda-c);
    } else {
        lambda = zeros(be.n,1);
        return -ginvH*c;
    }
}

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
 * 调用方法：sol=quaprog(H,c,Ae,be,Ai,bi,x0,[MAXN])
 * ****************************************************************/
Matrix quaprog(const Matrix &H, const Matrix &c, const Matrix &Ae, const Matrix &be, 
               const Matrix &Ai, const Matrix &bi, const Matrix &x0, const int MAXN=1000, const double eps=1e-10){
    Matrix x = x0;
    int n = x.n, ne = be.n, ni = bi.n, step = 0;
    bool *idx = new bool[ni];
    for(int i = 0; i < ni; i++)
        idx[i] = !( value(Ai.getSubmatrix(i,i,0,n-1)*x) > bi[i][0]+eps );
    for(; step <= MAXN; step++){
        Matrix Aee;
        if(ne) Aee = Ae;
        for(int j = 0; j < ni; j++)
            if(idx[j]) Aee = mergeRow(Aee,Ai.getSubmatrix(j,j,0,n-1));
        Matrix gk = H*x+c, lamk, dk;
        int m = Aee.n;
        dk = quaprog_subproblem(H,gk,Aee,zeros(m,1),lamk);
        if(dk.vecnorm(2) <= eps){
            double y = 0.0; int jk;
            if(lamk.n > ne)
                for(int j = ne; j < lamk.n; j++)
                    if(lamk[j][0]<y) y=lamk[j][0], jk=j;
            if(y >= 0) break;
            for(int i = 0; i < ni; i++)
                if(idx[i] && ne+checkednum(idx,0,i)==jk+1){
                    idx[i] = false;
                    break;
                }
        } else {
            int ti;
            double alpha = 1.0, tm = 1.0;
            for(int i = 0; i < ni; i++){
                Matrix tmp = Ai.getSubmatrix(i,i,0,n-1);
                if( !idx[i] && value(tmp*dk)<0 ){
                    double tm1 = (bi[i][0]-value(tmp*x))/value(tmp*dk);
                    if(tm1 < tm) tm = tm1, ti = i;
                }
            }
            alpha = std::min(alpha, tm);
            x = x+alpha*dk;
            if(tm < 1) idx[ti] = true;
        }
    }
#ifndef SILENCE
    std::cerr << "---------- Quadric Programming Active Set Method ----------" << std::endl;
    if(step<=MAXN) std::cerr << "Finished Succesfully. Total Steps: " << step << std::endl;
    else std::cerr << "Terminated. Too many steps." << std::endl;
    std::cerr << "Optimal Point: " << x.T() << std::endl;
    std::cerr << "OPtimal Value: " << 0.5*value(x.T()*H*x)+value(c.T()*x) << std::endl << std::endl;
#endif
    return x;
}

#endif