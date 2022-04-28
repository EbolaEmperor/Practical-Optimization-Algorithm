/*****************************************************************************************
 *
 * 这是一个BFGS拟牛顿法的通用最优化程序，可自动求导，依赖于：matrix.h
 * 用法：sol=bfgs(f,x0,[err],[rho])
 * 其中f是待求解函数，x0是初始迭代位置，err是允许误差，rho是AG条件中的参数
 * 返回值是BFGS方法求得的最小值点
 * f和grad的定义方法参照本例程
 * 
 * 注：本程序的一维搜索采用Goldstein准则
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 ****************************************************************************************/

#include "../matrix.h"
#include <cmath>
using namespace std;

// 以下是BFGS拟牛顿法通用子程序部分，不必修改
//********************************************************************************

double armijo_goldstein(double (*f)(const Matrix&), const Matrix &initial, 
                        const Matrix &direction, const double rho=0.4, const double t=100){
    double a = 0, b, alpha = 1;
    double f0 = f(initial), g0=value(direction.T()*gradient(f,initial));
    bool b_is_inf = true;
    int maxstep = 200;
    while(1){
        double fa = f(initial + alpha*direction);
        //cout << a << " " << b << " " << alpha << " " << fa << " " << f0+rho*alpha*g0 << " " << f0+(1-rho)*alpha*g0 << endl;
        if(!--maxstep) return alpha;
        if(fa <= f0+rho*alpha*g0){
            if(fa >= f0+(1-rho)*alpha*g0) break;
            else{
                if(b_is_inf) alpha = t*alpha;
                else alpha = (a+b)/2;
            }
        } else {
            b = alpha;
            b_is_inf = false;
            alpha = (a+b)/2;
        }
    }
    return alpha;
}

Matrix bfgs(double (*f)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.4){
    const int n = current.n;
    Matrix B = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = armijo_goldstein(f,current,direction,rho);
        Matrix s = alpha*direction;
        Matrix y = gradient(f,current+s) - gradient(f,current);
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

int main(){
    cin >> n;
    Matrix x(n,1);
    for(int i = 0; i < n; i++)
        x[i][0] = (i&1) ? 1 : -1.2;
    x = bfgs(f, x, 1e-5, 0.1);
    cout << "min f = f(" << x.T() << ") = " << f(x) << endl;
    return 0;
}