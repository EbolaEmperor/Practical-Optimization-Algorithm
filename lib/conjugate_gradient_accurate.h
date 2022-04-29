#ifndef _CONGUGATE_GRADIENT_ACCURATE_
#define _CONGUGATE_GRADIENT_ACCURATE_

#include <iostream>
#include "matrix_fraction.h"

/***************************************************************
 *
 * 这是一个用朴素共轭梯度法精确求解方程Ax=b的通用最优化程序
 * 调用CG_accurate(A,b,x)即可，其中x为初始迭代位置
 * 该算法能在有限步收敛到精确解，但速度很慢
 *
 **************************************************************/
fracMatrix CG_accurate(const fracMatrix &A, const fracMatrix &b, fracMatrix x){
    fracMatrix r = A*x-b;
    fracMatrix p = -r;
    long long step = 0;
    while(!r.iszero()){
        step++;
        std::cout << "Step=" << step << "   x=" << x.T() << std::endl;
        fraction alpha = - value(r.T()*p) / value(p.T()*A*p);
        x = x + alpha*p;
        fraction tmp = value(r.T()*r);
        r = r + alpha*A*p;
        fraction beta = value(r.T()*r) / tmp;
        p = -r + beta*p;
    }
    std::cout << "Steps: " << step << std::endl;
    return x;
}

#endif