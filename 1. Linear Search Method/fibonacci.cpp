#include <bits/stdc++.h>
using namespace std;

double argmin_fibonacci(double (*f)(double), double a, double b, double eps=1e-6)
//  计算一维函数极小值点的Fibonacci分割法
{
    double f0=1, f1=2, f2=3;
    while(b-a > eps)
    {
        double m1 = a+f0/f2*(b-a);
        double m2 = a+f1/f2*(b-a);
        if(f(m1)<f(m2)) b=m2;
        else a=m1;
        f0=f1;
        f1=f2;
        f2=f0+f1;
    }
    return a;
}

double f(double x)
{
    return x+1/x+x*x+0.5*exp(x);
}

int main()
{
    double a = argmin_fibonacci(f, 0, 2);
    printf("min f = f(%lf) = %lf\n", a, f(a));
    return 0;
}
