/********************************************************************************************************
 *
 * 这是一个DFP拟牛顿法的通用最优化程序，可自动求导，依赖于：matrix.h
 * 用法：sol=dfp(f,x0,[err],[rho],[sigma])
 * 其中f是待求解函数，x0是初始迭代位置，err是允许误差，rho,sigma是简单准则中的参数
 * 返回值是DFP方法求得的最小值点
 * f和grad的定义方法参照本例程
 * 
 * 注：本程序的一维搜索采用Wolfe准则
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 ********************************************************************************************************/

#include "../matrix.h"
#include <cmath>
using namespace std;

// 以下是BFGS拟牛顿法通用子程序部分，不必修改
//********************************************************************************

double zoom(double (*f)(const Matrix&), const Matrix &initial, const Matrix & direction, 
            double a, double b, const double rho, const double sigma){
    double alpha, fa, ga;
    double f0 = f(initial), g0 = value(direction.T()*gradient(f,initial));
    int step = 0;
    while(b-a>1e-8 && ++step<100){
        alpha = (a+b)/2;
        fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = value(direction.T()*gradient(f,initial+alpha*direction));
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell(double (*f)(const Matrix&), const Matrix &initial, const Matrix & direction, 
                    const double rho=0.4, const double sigma=0.7){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(initial), g0 = value(direction.T()*gradient(f,initial));
    int step = 0;
    while(++step<100){
        double fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            return zoom(f, initial, direction, a, alpha, rho, sigma);
        double ga = value(direction.T()*gradient(f,initial+alpha*direction));
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom(f, initial, direction, a, alpha, rho, sigma);
        a = alpha;
        alpha = (a+b)/2;
    }
    return alpha;
}

Matrix dfp(double (*f)(const Matrix&), Matrix current, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    const int n = current.n;
    Matrix H = eye(n);
    int step = 0;
    while(gradient(f,current).vecnorm(2) > err){
        step++;
        Matrix direction = -(H*gradient(f,current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = gradient(f,current+s) - gradient(f,current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        H = H - (1.0/value(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/value(s.T()*s)) * (s*s.T());
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
    x = dfp(f, x, 1e-5, 0.1, 0.4);
    cout << "min f = f(" << x.T() << ") = " << f(x) << endl;
    return 0;
}