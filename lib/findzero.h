#pragma once

#include <complex>
#include <iostream>

using Complex = std::complex<double>;

template <typename Func, typename Num>
Num findzero_newton(Func f, Func df, Num x0, double tol = 1e-6, int max_iter = 1000) {
    Num x = x0;
    for (int i = 0; i < max_iter; ++i) {
        Num fx = f(x), dfx = df(x);
        if (dfx == 0.0)
            throw std::runtime_error("Derivative is zero, cannot proceed.");
        Num x_new = x - fx / dfx;
#ifdef DEBUG
        std::cout << "Iteration " << i << ": x = " << x 
                  << ", f(x) = " << fx << ", df(x) = " << dfx
                  << ", diff = " << std::abs(x_new - x) << std::endl;
#endif
        if (std::abs(x_new - x) < tol) return x_new;
        x = x_new;
    }
    std::cerr << "Warning: Maximum iterations reached without convergence." << std::endl;
    return x;
}

template <typename Func, typename Num>
Num findzero_newton_gradfree(Func f, Num x0, double tol = 1e-6, int max_iter = 1000) {
    Num x = x0;
    auto df = [f, tol](Num x) {
        const Num h = tol / 10; // Small step for numerical derivative
        return (f(x + h) - f(x - h)) / (2.0 * h);
    };
    for (int i = 0; i < max_iter; ++i) {
        Num fx = f(x), dfx = df(x);
        if (dfx == 0.0)
            throw std::runtime_error("Derivative is zero, cannot proceed.");
        Num x_new = x - fx / dfx;
#ifdef DEBUG
        std::cout << "Iteration " << i << ": x = " << x 
                  << ", f(x) = " << fx << ", df(x) = " << dfx
                  << ", diff = " << std::abs(x_new - x) << std::endl;
#endif
        if (std::abs(x_new - x) < tol) return x_new;
        x = x_new;
    }
    std::cerr << "Warning: Maximum iterations reached without convergence." << std::endl;
    return x;
}

template <typename Func, typename Num>
Num findzero_secant(Func f, Num x0, Num x1,
                    double tol = 1e-6, int max_iter = 1000) {
    Num f0 = f(x0), f1 = f(x1);
    for (int i = 0; i < max_iter; ++i) {
        if (f1 - f0 == 0.0)
            throw std::runtime_error("Function values are too close, cannot proceed.");
        Num x_new = x1 - f1 * (x1 - x0) / (f1 - f0);
        Num f_new = f(x_new);
#ifdef DEBUG
        std::cout << "Iteration " << i << ": x = " << x_new
                  << ", f(x) = " << f_new << ", diff = " << std::abs(x_new - x1) << std::endl;
#endif
        if (std::abs(f_new) < tol) return x_new;
        x0 = x1; f0 = f1;
        x1 = x_new; f1 = f_new;
    }
    std::cerr << "Warning: Maximum iterations reached without convergence." << std::endl;
    return x1;
}

double _sign(double x) {
    if (fabs(x) < 1e-15) return 0.0;
    else return x > 0.0 ? 1.0 : -1.0;
}

template <typename Func>
Complex findzero_muller(Func f, Complex x0, Complex x1, Complex x2,
                        double tol = 1e-6, int max_iter = 1000) {
    for (int i = 0; i < max_iter; ++i) {
        Complex f0 = f(x0), f1 = f(x1), f2 = f(x2);
        Complex df12 = (f2 - f1) / (x2 - x1);
        Complex df01 = (f1 - f0) / (x1 - x0);
        Complex a = (df12 - df01) / (x2 - x0);
        Complex b = df12 + a * (x2 - x1);
        Complex c = f2;
        Complex discriminant = b * b - 4.0 * a * c;
        Complex x3 = x2 - 2.0 * c * _sign(b.real()) / (std::sqrt(discriminant) + std::abs(b));
        Complex f3 = f(x3);
#ifdef DEBUG
        std::cout << "Iteration " << i << ": x = " << x3
                  << ", f(x) = " << f3 << ", diff = " << std::abs(x3 - x2) << std::endl;
        std::cout << "             a = " << a << ", b = " << b << ", c = " << c
                  << ", discriminant = " << discriminant << std::endl;
#endif
        if (std::abs(f3) < tol) return x3;
        x0 = x1; x1 = x2; x2 = x3;
    }
    return x2;
}