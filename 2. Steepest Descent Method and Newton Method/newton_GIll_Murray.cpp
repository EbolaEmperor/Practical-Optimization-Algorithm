/***************************************************************
 *
 * 这是一个Gill-Murray修正牛顿法的通用最优化程序
 * 其中线搜索采用了牛顿单点插值法
 * 实际应用时需要将：
 *   - f设为待优化函数
 *   - grad设为f的梯度
 *   - hessian设为f的Hessian矩阵
 *   - current设为初始迭代位置
 * 其它模块无需作修改，依赖于matrix.h
 * 例程中的函数是Beale函数
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/


#include <iostream>
#include <cstdio>
#include "../matrix.h"
using namespace std;

// 用户自定义部分
//*********************************************************************

double f(const Matrix & x)
// 返回f在x处的函数值
{
    return (1.5-x[0][0]+x[0][0]*x[1][0])*(1.5-x[0][0]+x[0][0]*x[1][0]) + 
           (2.25-x[0][0]+x[0][0]*x[1][0]*x[1][0])*(2.25-x[0][0]+x[0][0]*x[1][0]*x[1][0]) +
           (2.625-x[0][0]+x[0][0]*x[1][0]*x[1][0]*x[1][0])*(2.625-x[0][0]+x[0][0]*x[1][0]*x[1][0]*x[1][0]);
}

Matrix grad(const Matrix & x)
// 返回f在x处的梯度向量
{
    Matrix g(2,1);
    g[0][0] = 2*(x[1][0]-1)*(1.5-x[0][0]+x[0][0]*x[1][0]) + 
              2*(x[1][0]*x[1][0]-1)*(2.25-x[0][0]+x[0][0]*x[1][0]*x[1][0]) +
              2*(x[1][0]*x[1][0]*x[1][0]-1)*(2.625-x[0][0]+x[0][0]*x[1][0]*x[1][0]*x[1][0]);
    g[1][0] = 2*x[0][0]*(1.5-x[0][0]+x[0][0]*x[1][0]) + 
              4*x[0][0]*x[1][0]*(2.25-x[0][0]+x[0][0]*x[1][0]*x[1][0]) +
              6*x[0][0]*x[1][0]*x[1][0]*(2.625-x[0][0]+x[0][0]*x[1][0]*x[1][0]*x[1][0]);
    return g;
}

Matrix hessian(const Matrix &x)
// 返回f在x处的Hessian矩阵
{
    Matrix h(2,2);
    h[0][0] = (2*(x[1][0]-1)*(x[1][0]-1) + 
              2*(x[1][0]*x[1][0]-1)*(x[1][0]*x[1][0]-1) +
              2*(x[1][0]*x[1][0]*x[1][0]-1)*(x[1][0]*x[1][0]*x[1][0]-1));
    h[0][1] = h[1][0] = (2*(1.5-x[0][0]+x[0][0]*x[1][0]) + 2*x[0][0]*(x[1][0]-1) +
                        4*x[1][0]*(2.25-x[0][0]+x[0][0]*x[1][0]*x[1][0]) + 4*x[0][0]*x[1][0]*(x[1][0]*x[1][0]-1) +
                        6*x[1][0]*x[1][0]*(2.625-x[0][0]+x[0][0]*x[1][0]*x[1][0]*x[1][0]) + 6*x[0][0]*x[1][0]*x[1][0]*(x[1][0]*x[1][0]*x[1][0]-1));
    h[1][1] = (2*x[0][0]*x[0][0] + 
              4*x[0][0]*(2.25-x[0][0]+x[0][0]*x[1][0]*x[1][0]) + 8*x[0][0]*x[0][0]*x[1][0]*x[1][0] + 
              12*x[0][0]*x[1][0]*(2.625-x[0][0]+x[0][0]*x[1][0]*x[1][0]*x[1][0]) + 18*x[0][0]*x[0][0]*x[1][0]*x[1][0]*x[1][0]*x[1][0]);
    return h;
}

// 通用程序部分
//*********************************************************************

// 返回f在x处沿方向direct的方向导数，其中direct是单位方向向量
double directionalDerivative(const Matrix & x, const Matrix & direct){
    Matrix g = grad(x);
    return value(g.T()*direct);
}

// 由最速下降法确定的一维搜索方向，是单位方向向量
Matrix searchDirection;

// 最速下降法当前迭代到的位置
Matrix current;

// 返回f(a+td)，其中a是最速下降法当前迭代到的位置，d是由最速下降法确定的一维搜索方向
double f_1dim(double t){
    return f(current+t*searchDirection);
}

// 返回f(a+td)关于t的导数
double fder_1dim(double t){
    return directionalDerivative(current+t*searchDirection, searchDirection);
}

// 返回f(a+td)关于t的二阶导数
double fder2_1dim(double t){
    return value( searchDirection.T() * hessian(current+t*searchDirection) * searchDirection) ;
}

// 计算一维函数极小值点的Newton单点二次插值法，需要传入导函数、二阶导函数
double argmin_1dim(double (*f1)(double), double (*f2)(double), double x, double eps=1e-6){
    int step = 0;
    while(1){
        double y = x-f1(x)/f2(x);
        if(fabs(y-x)<eps || ++step >= 100) return y;
        // 注：这里使用了一个小小的Trick，限制线搜索的迭代步数最多100步，以提高效率
        else x=y;
    }
}

// Gill-Murray修正牛顿法的主程序
double min_newton_GM(double eps = 1e-6){
    int step = 0;
    while(grad(current).vecnorm(2)>eps){
        step++;
        searchDirection = solveByLDL_GM(hessian(current), -grad(current));
        searchDirection = (1.0/searchDirection.vecnorm(2)) * searchDirection;
        double lambda = argmin_1dim(fder_1dim, fder2_1dim, 0, 1e-2*eps);
        current = current + lambda * searchDirection;
        cout << "step: " << step << "    current=" << current.T() << "    direction=" << searchDirection.T() << "    f=" << f(current) << endl;
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    cout << "Total steps: " << step << endl;
    return f(current);
}

// 主程序部分，在主程序中需要给出初始迭代向量
//*********************************************************************

int main(){
    current = Matrix(2,1);
    current[0][0] = -1.2;
    current[1][0] = 1.0;
    double ans = min_newton_GM(1e-3);
    printf("min f = f(%.3lf,%.3lf) = %.3lf\n", current[0][0], current[1][0], ans);
    return 0;
}