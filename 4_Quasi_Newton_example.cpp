#include "lib/quassi_newton.h"
#include "func/rosenbrock.h"
#include <cmath>
using namespace std;
using namespace rosenbrock;

// 以扩展Rosenbrock函数为例
//********************************************************************************

int main(){
    cin >> n;
    Matrix x=initial(), y;

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