#include <bits/stdc++.h>
using namespace std;

double argmin_newton(double (*f1)(double), double (*f2)(double), double x, double eps=1e-6)
//  计算一维函数极小值点的Newton单点二次插值法，需要传入导函数、二阶导函数
{
    while(1)
    {
        double y = x-f1(x)/f2(x);
        if(fabs(y-x)<eps) return y;
        else x=y;
    }
}

double f(double x)
{
    return x+1/x+x*x+0.5*exp(x);
}
double f1(double x)
{
    return 1-1/(x*x)+2*x+0.5*exp(x);
}
double f2(double x)
{
    return 2/(x*x*x)+2+0.5*exp(x);
}

int main()
{
    double a = argmin_newton(f1, f2, 2);
    printf("min f = f(%lf) = %lf\n", a, f(a));
    return 0;
}