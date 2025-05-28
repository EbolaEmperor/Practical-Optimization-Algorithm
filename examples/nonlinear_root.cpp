#define DEBUG
#include "lib/findzero.h"
using namespace std;

double f(double x) {
    return x * x * x + 10.0 * x - 20.0;
}

double df(double x) {
    return 3.0 * x * x + 10.0;
}

Complex compf(Complex x) {
    return x * x * x + 10.0 * x - 20.0;
}

int main(){
    try {
        cout << "------------------- Newton Method -------------------" << endl;
        double root = findzero_newton(f, df, 1.0, 1e-12);
        cout << "\nRoot found: " << root << endl;
        cout << "f(root) = " << f(root) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------- Newton Method (Grad-Free) -------------" << endl;
        double root = findzero_newton_gradfree(f, 1.0, 1e-12);
        cout << "\nRoot found: " << root << endl;
        cout << "f(root) = " << f(root) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------------- Secant Method -------------------" << endl;
        double root = findzero_secant(f, 1.0, 1.2, 1e-12);
        cout << "\nRoot found: " << root << endl;
        cout << "f(root) = " << f(root) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    try {
        cout << "------------------- Muller Method ------------------" << endl;
        Complex root = findzero_muller(compf, 
                                       Complex(1.5, 0.0), 
                                       Complex(1.75, 0.0), 
                                       Complex(2.0, 0.0), 
                                       1e-12);
        cout << "\nRoot found: " << root.real() << "+" << root.imag() << "i" << endl;
        cout << "|f(root)| = " << std::abs(compf(root)) << endl << endl;
    } catch (const std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}