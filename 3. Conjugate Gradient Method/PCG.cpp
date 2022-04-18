/***************************************************************
 *
 * 这是一个用预优共轭梯度法求解方程Ax=b的通用最优化程序
 * 调用PCG(A,b,x,M,err)即可，其中x为初始迭代位置, M是预优因子
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#include <bits/stdc++.h>
#include "../matrix.h"
using namespace std;

Matrix PCG(const Matrix &A, const Matrix &b, Matrix x, Matrix M, const double err = 1e-6){
    Matrix r = A*x-b, y = solve(M,r), p = -y;
    long long step = 0;
    while(r.vecnorm(2) >= err){
        step++;
        double alpha = value(r.T()*y) / value(p.T()*A*p);
        x = x + alpha*p;
        double tmp = value(r.T()*y);
        r = r + alpha*A*p;
        y = solve(M,r);
        double beta = value(r.T()*y) / tmp;
        p = -y + beta*p;
    }
    cout << "Steps: " << step << endl;
    return x;
}

int main(){

    return 0;
}