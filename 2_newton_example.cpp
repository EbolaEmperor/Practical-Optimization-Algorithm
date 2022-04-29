#include <iostream>
#include <cstdio>
#include "lib/newton.h"
#include "func/rosenbrock.h"
#include "func/beale.h"
using namespace std;

int main(){
    rosenbrock::n = 2;
    Matrix y;

    cout << "Pure Newton Method" << endl;
    y = newton(rosenbrock::f, rosenbrock::grad, rosenbrock::hessian, rosenbrock::initial(), 1e-6);
    cout << "min f = (" << y.T() << ") = " << rosenbrock::f(y) << endl << endl;


    cout << "Newton Method With Wolfe Search" << endl;
    y = newton_wolfe(rosenbrock::f, rosenbrock::grad, rosenbrock::hessian, rosenbrock::initial(), 1e-6, 0.1, 0.3);
    cout << "min f = (" << y.T() << ") = " << rosenbrock::f(y) << endl << endl;

    cout << "Gill-Murray Corrected Newton Method" << endl;
    y = newton_gillmurray(beale::f, beale::grad, beale::hessian, beale::initial(), 1e-6, 0.2, 0.3);
    // 注：这个例子中，迭代无法收敛
    cout << "min f = (" << y.T() << ") = " << beale::f(y) << endl << endl;
    return 0;
}