/***************************************************************
 *
 * 这是一个最速下降法的通用最优化程序
 * 其中线搜索采用了牛顿单点插值法
 * 实际应用时需要将：
 *   - f设为待优化函数
 *   - grad设为f的梯度
 *   - hessian设为f的Hessian矩阵
 *   - current设为初始迭代位置
 * 其它模块无需作修改
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/


#include <bits/stdc++.h>
using namespace std;

// 用户自定义部分
//*********************************************************************

double f(const vector<double> & x)
// 返回f在x处的函数值
{
    return 100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]) + (1-x[0])*(1-x[0]);
}

vector<double> grad(const vector<double> & x)
// 返回f在x处的梯度向量
{
    vector<double> g(2);
    g[0] = 400*x[0]*(x[0]*x[0]-x[1]) + 2*(x[0]-1);
    g[1] = 200*(x[1]-x[0]*x[0]);
    return g;
}

vector< vector<double> > hessian(const vector<double> &x)
// 返回f在x处的Hessian矩阵
{
    vector< vector<double> > h(2);
    h[0].push_back(1200*x[0]*x[0]+2);
    h[1].push_back(-400*x[0]);
    h[1].push_back(200);
    return h;
}

// 通用程序部分
//*********************************************************************

vector<double> operator + (const vector<double> &a, const vector<double> &b)
// 向量加法运算
{
    assert(a.size() == b.size());
    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++)
        c[i] = a[i] + b[i];
    return c;
}

vector<double> operator - (const vector<double> &b)
// 向量取负运算
{
    vector<double> c(b.size());
    for(int i = 0; i < b.size(); i++)
        c[i] = -b[i];
    return c;
}

vector<double> operator * (const double &k, const vector<double> &b)
// 向量数乘运算
{
    vector<double> c(b.size());
    for(int i = 0; i < b.size(); i++)
        c[i] = k * b[i];
    return c;
}

double operator * (const vector<double> &a, const vector<double> &b)
// 向量内积运算
{
    assert(a.size() == b.size());
    double c = 0;
    for(int i = 0; i < a.size(); i++)
        c += a[i] * b[i];
    return c;
}

double norm(const vector<double> &a)
// 向量的范数（二范数）
{
    double x = 0;
    for(double c : a)
        x += c*c;
    return sqrt(x);
}

double innerProduct(const vector< vector<double> > & A, const vector<double> &x)
// 计算x^TAx的值，其中A是一个Hessian阵
{
    assert(x.size() == A.size());
    double res = 0;
    for(int i = 0; i < x.size(); i++)
    {
        res += A[i][i]*x[i]*x[i];
        for(int j = 0; j < i-1; j++)
            res += 2*A[i][j]*x[i]*x[j];
    }
    return res;
}

double directionalDerivative(const vector<double> & x, const vector<double> & direct)
// 返回f在x处沿方向direct的方向导数，其中direct是单位方向向量
{
    vector<double> g = grad(x);
    return g*direct;
}

vector<double> searchDirection;
// 由最速下降法确定的一维搜索方向，是单位方向向量

vector<double> current;
// 最速下降法当前迭代到的位置

double f_1dim(double t)
// 返回f(a+td)，其中a是最速下降法当前迭代到的位置，d是由最速下降法确定的一维搜索方向
{
    return f(current+t*searchDirection);
}

double fder_1dim(double t)
// 返回f(a+td)关于t的导数
{
    return directionalDerivative(current+t*searchDirection, searchDirection);
}

double fder2_1dim(double t)
// 返回f(a+td)关于t的二阶导数
{
    return innerProduct(hessian(current+t*searchDirection), searchDirection);
}

double argmin_1dim(double (*f1)(double), double (*f2)(double), double x, double eps=1e-6)
// 计算一维函数极小值点的Newton单点二次插值法，需要传入导函数、二阶导函数
{
    while(1)
    {
        double y = x-f1(x)/f2(x);
        if(fabs(y-x)<eps) return y;
        else x=y;
    }
}

double min_SDM(double eps = 1e-6)
// 最速下降法的主程序
{
    int step = 0;
    while(norm(grad(current))>eps)
    {
        searchDirection = -grad(current);
        searchDirection = (1.0/norm(searchDirection)) * searchDirection;
        double lambda = argmin_1dim(fder_1dim, fder2_1dim, 0, eps);
        current = current + lambda * searchDirection;
        printf("step: %d,  current: (%lf, %lf),  direction: (%lf, %lf),  f = %lf\n", ++step, current[0], current[1], searchDirection[0], searchDirection[1], f(current));
    }
    return f(current);
}

int main()
{
    current.push_back(-1.2);
    current.push_back(1.0);
    double ans = min_SDM(1e-3);
    printf("min f = f(%.3lf,%.3lf) = %.3lf\n", current[0], current[1], ans);
    return 0;
}