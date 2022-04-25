/********************************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，依赖于：matrix.h
 * 用法：sol=bfgs(f,grad,n,x0,[err],[rho],[sigma])
 * 其中f是待求解函数，grad是待求解函数的梯度，n是函数的维数，x0是初始迭代位置，err是允许误差，rho,sigma是简单准则中的参数
 * 返回值是BFGS方法求得的最小值点
 * f和grad的定义方法参照本例程
 * 
 * 注：本程序的一维搜索采用简单准则
 * 简单准则是指满足：f(xk + rho^m * dk) <= f(xk) + sigma * rho^m * gk^T * dk 的最小m
 * 找到m后令 x_k+1 = xk + rho^m *dk 进入下一次迭代
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 ********************************************************************************************************/

#include "../matrix.h"
#include <cmath>
using namespace std;

// 以下是BFGS拟牛顿法通用子程序部分，不必修改
//********************************************************************************

double simple_search(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const Matrix &initial, const Matrix & direction, 
                    const double rho, const double sigma){
    int step = 0;
    double alpha = 1, f0 = f(initial), g0 = value(direction.T()*grad(initial));
    while(++step < 50){
        double fa = f(initial + alpha * direction);
        if(fa <= f0 + sigma * alpha * g0)
            return alpha;
        alpha = alpha * rho;
    }
    return alpha;
}

Matrix bfgs(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const int &n, Matrix current, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    Matrix B = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = simple_search(f,grad,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
    }
    cout << "Total Steps: " << step << endl;
    return current;
}

// 以下是自定义区域，以扩展Rosenbrock函数为例
//********************************************************************************

int n;

double f(const Matrix &x){
    double res = 0;
    for(int i = 0; i < n-1; i++)
        res += 100*(x[i+1][0]-x[i][0]*x[i][0])*(x[i+1][0]-x[i][0]*x[i][0]) + (1-x[i][0])*(1-x[i][0]);
    return res;
}

Matrix grad(const Matrix &x){
    Matrix g(n,1);
    for(int i = 0; i < n-1; i++){
        g[i][0]   += 400*x[i][0]*(x[i][0]*x[i][0]-x[i+1][0]) + 2*(x[i][0]-1);
        g[i+1][0] += 200*(x[i+1][0]-x[i][0]*x[i][0]);
    }
    return g;
}

int main(){
    cin >> n;
    Matrix x(n,1);
    for(int i = 0; i < n; i++)
        x[i][0] = (i&1) ? 1 : -1.2;
    x = bfgs(f, grad, n, x, 1e-5, 0.1, 0.1);
    cout << "min f = f(" << x.T() << ") = " << f(x) << endl;
    return 0;
}