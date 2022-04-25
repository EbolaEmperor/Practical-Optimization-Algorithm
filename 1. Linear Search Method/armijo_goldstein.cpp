#include <bits/stdc++.h>
using namespace std;

double armijo_goldstein(double (*f)(double), double (*g)(double), const double rho=0.4, const double t=100){
//  Armijo-Goldstein不精确一维搜索方法；需要传入函数、导函数；Armijo-Goldstein准则中的参数rho可选，默认0.4；步进因子t可选，默认100
    double a = 0, b, alpha = 1;
    double f0 = f(0), g0 = g(0);
    bool b_is_inf = true;
    while(1){
        double fa = f(alpha);
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

double f(double x)
{
    x += 0.2;
    return x+1/x+x*x+0.5*exp(x);
}
double f1(double x)
{
    x += 0.2;
    return 1-1/(x*x)+2*x+0.5*exp(x);
}

int main()
{
    double a = armijo_goldstein(f, f1);
    printf("accepted min f = f(%lf) = %lf\n", a, f(a));
    return 0;
}