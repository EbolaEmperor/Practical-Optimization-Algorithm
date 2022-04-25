#include <bits/stdc++.h>
using namespace std;

double zoom(double (*f)(double), double (*g)(double), double a, double b, const double rho, const double sigma){
    double alpha, fa, ga;
    double f0 = f(0), g0 = g(0);
    int step = 0;
    while(b-a>1e-8 && ++step<100){
        alpha = (a+b)/2;
        fa = f(alpha);
        if(fa > f0+rho*alpha*g0)
            b = alpha;
        else{
            ga = g(alpha);
            if(fabs(ga) <= -sigma*g0) return alpha;
            else if(ga*(b-a) >= 0) b = alpha;
            else a = alpha;
        }
    }
    return alpha;
}

double wolfe_powell(double (*f)(double), double (*g)(double), const double rho=0.4, const double sigma=0.7){
//  Wolfe-Powell不精确一维搜索方法；需要传入函数、导函数；Wolfe-Powell准则中的参数rho, sigma可选，默认0.4, 0.7；
    double a = 0, b = 1e3, alpha = 1;
    double f0 = f(0), g0 = g(0);
    int step = 0;
    while(++step<100){
        double fa = f(alpha);
        if(fa > f0+rho*alpha*g0)
            return zoom(f, g, a, alpha, rho, sigma);
        double ga = g(alpha);
        if(fabs(ga) <= -sigma*g0)
            return alpha;
        if(ga >= 0)
            return zoom(f, g, a, alpha, rho, sigma);
        a = alpha;
        alpha = (a+b)/2;
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
    double a = wolfe_powell(f, f1);
    printf("accepted min f = f(%lf) = %lf\n", a, f(a));
    return 0;
}