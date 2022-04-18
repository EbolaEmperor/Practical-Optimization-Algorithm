/***************************************************************
 *
 * 这是一个用朴素共轭梯度法精确求解方程Ax=b的通用最优化程序
 * 调用CG_accurate(A,b,x)即可，其中x为初始迭代位置
 * main函数中是使用示例，以n阶Hilbert矩阵为例
 * 该算法能在有限步收敛到精确解，但速度很慢
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#include <bits/stdc++.h>
#include "../matrix_fraction.h"
using namespace std;

fracMatrix CG_accurate(const fracMatrix &A, const fracMatrix &b, fracMatrix x){
    fracMatrix r = A*x-b;
    fracMatrix p = -r;
    long long step = 0;
    while(!r.iszero()){
        step++;
        cout << "Step=" << step << "   x=" << x.T() << endl;
        fraction alpha = - value(r.T()*p) / value(p.T()*A*p);
        x = x + alpha*p;
        fraction tmp = value(r.T()*r);
        r = r + alpha*A*p;
        fraction beta = value(r.T()*r) / tmp;
        p = -r + beta*p;
    }
    cout << "Steps: " << step << endl;
    return x;
}

int main(){
    int n;
    cin >> n;
    fracMatrix x = CG_accurate(hilbert(n), ones(n,1), zeros(n,1));
    cout << "ans = " << x.T() << endl;
    return 0;
}
