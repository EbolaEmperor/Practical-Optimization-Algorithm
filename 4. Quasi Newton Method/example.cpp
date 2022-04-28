#include "quassi_newton.h"
#include <cmath>
using namespace std;

// 以扩展Rosenbrock函数为例
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
    Matrix x(n,1), y;
    for(int i = 0; i < n; i++)
        x[i][0] = (i&1) ? 1 : -1.2;

    cout << "BFGS Method (with Wolfe)" << endl;
    y = bfgs(f, grad, x, 1e-5, 0.1, 0.4);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl << endl;

    cout << "BFGS Method (with Goldstein)" << endl;
    y = bfgs_goldstein(f, grad, x, 1e-5, 0.2);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl << endl;

    cout << "BFGS Method (with Simple)" << endl;
    y = bfgs_simple(f, grad, x, 1e-5, 0.1, 0.2);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl << endl;

    cout << "DFP Method (with Wolfe)" << endl;
    y = dfp(f, grad, x, 1e-5, 0.1, 0.4);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl << endl;

    cout << "Broyden Method (with Wolfe)" << endl;
    y = broyden(f, grad, x, 0.8, 1e-5, 0.1, 0.4);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl << endl;
    return 0;
}