#include <iostream>
#include "lib/conjugate_gradient.h"
#include "func/rosenbrock.h"
using namespace std;

int main(){
    int n;
    cin >> n;

    cout << "CG Method solve Hilbert(n) Equation" << endl;
    Matrix x = CG(hilbert(n), ones(n,1), zeros(n,1), 1e-6);
    cout << "ans = " << x.T() << endl << endl;

    cout << "Restart FR-CG Method solve n-dim rosenbrock minimal" << endl;
    rosenbrock::n = n;
    Matrix x2 = CG_FR(rosenbrock::f, rosenbrock::grad, rosenbrock::initial(), 1e-6, 0.1, 0.4);
    cout << "min f = f(" << x2.T() << ") = " << rosenbrock::f(x2) << endl << endl;

    cout << "Restart FR-CG Method (gradfree) solve n-dim rosenbrock minimal" << endl;
    rosenbrock::n = n;
    Matrix x3 = CG_FR_gradfree(rosenbrock::f, rosenbrock::initial(), 1e-6, 0.1, 0.4);
    cout << "min f = f(" << x3.T() << ") = " << rosenbrock::f(x3) << endl << endl;
    return 0;
}