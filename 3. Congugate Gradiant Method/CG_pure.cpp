/***************************************************************
 *
 * 这是一个用朴素共轭梯度法求解方程Ax=b的通用最优化程序
 * 调用CG(A,b,x,err)即可，其中x为初始迭代位置
 * main函数中是使用示例，以n阶Hilbert矩阵为例
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#include <bits/stdc++.h>
#include "../matrix.h"
using namespace std;

Matrix CG(const Matrix &A, const Matrix &b, Matrix x, const double err = 1e-6){
    Matrix r = A*x-b;
    Matrix p = -r;
    long long step = 0;
    while(r.vecnorm(2) >= err){
        step++;
        double alpha = - value(r.T()*p) / value(p.T()*A*p);
        x = x + alpha*p;
        double tmp = value(r.T()*r);
        r = r + alpha*A*p;
        double beta = value(r.T()*r) / tmp;
        p = -r + beta*p;
    }
    cout << "Steps: " << step << endl;
    return x;
}

int main(){
    int n;
    cin >> n;
    Matrix x = CG(hilbert(n), ones(n,1), zeros(n,1), 1e-6);
    cout << "ans = " << x.T() << endl;
    return 0;
}