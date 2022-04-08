/***************************************************************
 *
 * 这是一个Gill-Murray修正牛顿法的通用最优化程序
 * 其中线搜索采用了牛顿单点插值法
 * 实际应用时需要将：
 *   - f设为待优化函数
 *   - grad设为f的梯度
 *   - hessian设为f的Hessian矩阵
 *   - current设为初始迭代位置
 * 其它模块无需作修改
 * 例程中的函数是Beale函数
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
    return (1.5-x[0]+x[0]*x[1])*(1.5-x[0]+x[0]*x[1]) + 
           (2.25-x[0]+x[0]*x[1]*x[1])*(2.25-x[0]+x[0]*x[1]*x[1]) +
           (2.625-x[0]+x[0]*x[1]*x[1]*x[1])*(2.625-x[0]+x[0]*x[1]*x[1]*x[1]);
}

vector<double> grad(const vector<double> & x)
// 返回f在x处的梯度向量
{
    vector<double> g(2);
    g[0] = 2*(x[1]-1)*(1.5-x[0]+x[0]*x[1]) + 
           2*(x[1]*x[1]-1)*(2.25-x[0]+x[0]*x[1]*x[1]) +
           2*(x[1]*x[1]*x[1]-1)*(2.625-x[0]+x[0]*x[1]*x[1]*x[1]);
    g[1] = 2*x[0]*(1.5-x[0]+x[0]*x[1]) + 
           4*x[0]*x[1]*(2.25-x[0]+x[0]*x[1]*x[1]) +
           6*x[0]*x[1]*x[1]*(2.625-x[0]+x[0]*x[1]*x[1]*x[1]);
    return g;
}

vector< vector<double> > hessian(const vector<double> &x)
// 返回f在x处的Hessian矩阵
{
    vector< vector<double> > h(2);
    h[0].push_back(2*(x[1]-1)*(x[1]-1) + 
                   2*(x[1]*x[1]-1)*(x[1]*x[1]-1) +
                   2*(x[1]*x[1]*x[1]-1)*(x[1]*x[1]*x[1]-1));
    h[1].push_back(2*(1.5-x[0]+x[0]*x[1]) + 2*x[0]*(x[1]-1) +
                   4*x[1]*(2.25-x[0]+x[0]*x[1]*x[1]) + 4*x[0]*x[1]*(x[1]*x[1]-1) +
                   6*x[1]*x[1]*(2.625-x[0]+x[0]*x[1]*x[1]*x[1]) + 6*x[0]*x[1]*x[1]*(x[1]*x[1]*x[1]-1));
    h[1].push_back(2*x[0]*x[0] + 
                   4*x[0]*(2.25-x[0]+x[0]*x[1]*x[1]) + 8*x[0]*x[0]*x[1]*x[1] + 
                   12*x[0]*x[1]*(2.625-x[0]+x[0]*x[1]*x[1]*x[1]) + 18*x[0]*x[0]*x[1]*x[1]*x[1]*x[1]);
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

vector< vector<double> > transpose(const vector< vector<double> > &A)
// 矩阵转置，仅能用于下三角阵转置为上三角阵
{
    int n = A.size();
    vector< vector<double> > B(n);
    for(int i = 0; i < n; i++)
        for(int j = i; j < n; j++)
            B[i].push_back(A[j][i]);
    return B;
}

vector< vector<double> > cholesky_GM(vector< vector<double> > A)
// Gill-Murray修正cholesky分解
{
    int n = A.size();
    double gamma=0, xi=0;
    for(int i = 0; i < n; i++)
    {
        gamma = max(gamma, fabs(A[i][i]));
        for(int j = 0; j < i; j++)
            xi = max(xi, fabs(A[i][j]));
        for(int j = i+1; j < n; j++)
            A[i].push_back(A[j][i]);
    }
    double nu = max(1.0, sqrt(n*n-1));
    double beta2 = max(max(gamma,xi/nu),1e-6);
    vector< vector<double> > c(n);
    for(int i = 0; i < n; i++)
        c[i].resize(n), c[i][i] = A[i][i];
    vector< vector<double> > L(n);
    for(int i = 0; i < n; i++)
        L[i].resize(i+1);
    for(int j = 0; j < n; j++)
    {
        int q = j;
        for(int k = j+1; k < n; k++)
            if(fabs(c[k][k])>fabs(c[q][q])) q = k;
        swap(A[j], A[q]);
        for(int k = 0; k < n; k++)
            swap(A[k][j], A[k][q]);
        for(int k = 0; k < j; k++)
            L[j][k] = c[j][k]/L[k][k];
        for(int i = j+1; i < n; i++)
        {
            c[i][j] = A[i][j];
            for(int k = 0; k < j; k++)
                c[i][j] -= c[i][k]*L[j][k];
        }
        double theta = 0;
        for(int k = j+1; k < n; k++)
            theta = max(theta, fabs(c[k][j]));
        L[j][j] = max(max(fabs(c[j][j]),theta*theta/beta2),1e-3);
        for(int i = j+1; i < n; i++)
            c[i][i] -= c[i][j]*c[i][j]/L[j][j];
    }
    return L;
}

vector<double> solveUnitLowerTriangular(const vector< vector<double> > &A, const vector<double> &b)
// 求解单位下三角线性方程组
{
    int n = A.size();
    vector<double> x;
    for(double v : b)
        x.push_back(v);
    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++)
            x[j] -= x[i] * A[j][i];
    return x;
}

vector<double> solveUnitUpperTriangular(const vector< vector<double> > &A, const vector<double> &b)
// 求解单位上三角线性方程组
{
    int n = A.size();
    vector<double> x;
    for(double v : b)
        x.push_back(v);
    for(int i = n-1; i >= 0; i--)
        for(int j = 0; j < i; j++)
            x[j] -= x[i] * A[j][i];
    return x;
}

vector<double> solveHessianEquation(const vector< vector<double> > &A, const vector<double> &b)
// 用Gill-Murray修改Cholesky分解法求解方程(A+E)x=b，其中A是一个Hessian阵（对称正定）
{
    vector< vector<double> > L = cholesky_GM(A);
    vector<double> y = solveUnitLowerTriangular(L, b);
    for(int i = 0; i < y.size(); i++)
        y[i] /= L[i][i];
    return solveUnitUpperTriangular(transpose(L), y);
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
    int step = 0;
    while(1)
    {
        double y = x-f1(x)/f2(x);
        if(fabs(y-x)<eps || ++step >= 100) return y;
        // 注：这里使用了一个小小的Trick，限制线搜索的迭代步数最多100步，以提高效率
        else x=y;
    }
}

double min_newton(double eps = 1e-6)
// 牛顿法的主程序
{
    int step = 0;
    while(norm(grad(current))>eps)
    {
        step++;
        searchDirection = solveHessianEquation(hessian(current), -grad(current));
        searchDirection = (1.0/norm(searchDirection)) * searchDirection;
        double lambda = argmin_1dim(fder_1dim, fder2_1dim, 0, 1e-2*eps);
        current = current + lambda * searchDirection;
        printf("step: %d,  current: (%lf, %lf),  direction: (%lf, %lf),  f = %lf\n", step, current[0], current[1], searchDirection[0], searchDirection[1], f(current));
        // 这是用于输出每一步迭代信息的测试语句，可以删除
    }
    printf("Total steps: %d\n", step);
    return f(current);
}

int main()
{
    current.push_back(-1.2);
    current.push_back(1.0);
    double ans = min_newton(1e-3);
    printf("min f = f(%.3lf,%.3lf) = %.3lf\n", current[0], current[1], ans);
    return 0;
}