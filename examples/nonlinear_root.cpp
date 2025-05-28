#define DEBUG
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
    return 0;
}