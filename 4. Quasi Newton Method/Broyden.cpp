/********************************************************************************************************
 *
 * 这是一个Broyden族校正拟牛顿法的通用最优化程序，依赖于：matrix.h
 * 用法：sol=broyden(f,grad,n,x0,[phi],[err],[rho],[sigma])
 * 其中f是待求解函数，grad是待求解函数的梯度，n是函数的维数，x0是初始迭代位置，err是允许误差，rho,sigma是Wolfe准则中的参数
 * phi是Broyden族校正中的重要参数，默认为1（即得BFGS校正），为0时就是DFP校正
 * 返回值是Broyden族校正方法求得的最小值点
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

double zoom(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const Matrix &initial, const Matrix & direction, 
            double a, double b, const double rho, const double sigma){
    double alpha, fa, ga;
    double f0 = f(initial), g0 = value(direction.T()*grad(initial));
    int step = 0;
    while(b-a>1e-8 && ++step<100){
        alpha = (a+b)/2;
        fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = value(direction.T()*grad(initial+alpha*direction));
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const Matrix &initial, const Matrix & direction, 
                    const double rho=0.4, const double sigma=0.7){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(initial), g0 = value(direction.T()*grad(initial));
    int step = 0;
    while(++step<100){
        double fa = f(initial + alpha*direction);
        if(fa > f0+rho*alpha*g0)
            return zoom(f, grad, initial, direction, a, alpha, rho, sigma);
        double ga = value(direction.T()*grad(initial+alpha*direction));
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom(f, grad, initial, direction, a, alpha, rho, sigma);
        a = alpha;
        alpha = (a+b)/2;
    }
    return alpha;
}

Matrix broyden(double (*f)(const Matrix&), Matrix (*grad)(const Matrix&), const int &n, Matrix current, 
               const double phi=1, const double err=1e-5, const double rho=0.55, const double sigma=0.4){
    Matrix B = eye(n), H = eye(n);
    int step = 0;
    while(grad(current).vecnorm(2) > err){
        step++;
        Matrix direction = -solveByLDL(B,phi*grad(current)) - (1-phi)*(H*grad(current));
        direction = (1.0/direction.vecnorm(2))*direction;
        double alpha = wolfe_powell(f,grad,current,direction,rho,sigma);
        Matrix s = alpha*direction;
        Matrix y = grad(current+s) - grad(current);
        current = current + s;
        //cout << "Step: " << step << "    current=" << current.T() << "    direction=" << direction.T() << "    f=" << f(current)  << "    ||grad||=" << grad(current).vecnorm(2) << endl;
        B = B - (1.0/value(s.T()*B*s)) * ((B*s)*(s.T()*B)) + (1.0/value(y.T()*s)) * (y*y.T());
        H = H - (1.0/value(y.T()*H*y)) * ((H*y)*(y.T()*H)) + (1.0/value(s.T()*y)) * (s*s.T());
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
    x = broyden(f, grad, n, x, 0.8, 1e-5, 0.1, 0.4);
    cout << "min f = f(" << x.T() << ") = " << f(x) << endl;
    return 0;
}