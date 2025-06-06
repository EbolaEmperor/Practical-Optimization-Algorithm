#include "lib/newton.h"
#include "lib/least_square.h"
#include "lib/findzero.h"
#include <iomanip>
using namespace std;

double f(double x) {
    return exp(x*x/2) - 3*sin(x);
}

double df(double x) {
    return x*exp(x*x/2) - 3*cos(x);
}

Complex compf(Complex x) {
    return x * x * x + 10.0 * x - 20.0;
}

ColVector g(const ColVector &x) {
    ColVector res(2);
    res[0] = x[0]*x[0] + x[1]*x[1] - 4;
    res[1] = x[0]*x[0]*x[0] - x[1];
    return res;
}

Matrix Jg(const ColVector &x) {
    Matrix res(2, 2);
    res[0][0] = 2 * x[0];
    res[0][1] = 2 * x[1];
    res[1][0] = 3 * x[0] * x[0];
    res[1][1] = -1;
    return res;
}

int main(){
    cout << std::setprecision(16);

    try {
        cout << "------------------- Newton Method -------------------" << endl;
        double root = findzero_newton(f, df, 2.0, 1e-13);
        cout << "\nRoot found: " << root << endl;
        cout << "f(root) = " << f(root) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------- Newton Method (Grad-Free) -------------" << endl;
        double root = findzero_newton_gradfree(f, 0.0, 1e-13);
        cout << "\nRoot found: " << root << endl;
        cout << "f(root) = " << f(root) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------------- Secant Method -------------------" << endl;
        double root = findzero_secant(f, 1.0, 1.2, 1e-13);
        cout << "\nRoot found: " << root << endl;
        cout << "f(root) = " << f(root) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    cout << std::setprecision(3);
    try {
        cout << "------------------- Muller Method ------------------" << endl;
        Complex root = findzero_muller(compf, 
                                       Complex(-1.0, 3.1), 
                                       Complex(-0.8, 3.3), 
                                       Complex(-0.6, 3.5), 
                                       1e-3);
        cout << "\nRoot found: " << root.real() << "+" << root.imag() << "i" << endl;
        cout << "|f(root)| = " << std::abs(compf(root)) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------------- Newton Method (2-dim) -------------------\n" << endl;
        ColVector x0(2);
        x0[0] = x0[1] = 1;
        ColVector root = newton_zero(g, Jg, x0, 1e-13);
        cout << "Root found: " << root.T() << endl;
        cout << "f(root) = " << g(root).T() << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------------- Newton Method Grad-Free (2-dim) -------------------\n" << endl;
        ColVector x0(2);
        x0[0] = x0[1] = 1;
        ColVector root = newton_zero_gradfree(g, x0, 1e-13);
        cout << "Root found: " << root.T() << endl;
        cout << "f(root) = " << g(root).T() << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    // Levenberg-Marquardt Method (2-dim)
    try {
        ColVector x0(2);
        x0[0] = x0[1] = 1;
        ColVector root = nonlinlsq_LM(g, Jg, x0, 1e-13);
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}