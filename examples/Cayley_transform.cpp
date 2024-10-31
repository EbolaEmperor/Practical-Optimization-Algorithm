// 寻找 Cayley transform 的例子，详见 Matrix Computation (Edition 4) 的习题 P5.1.8
#include "lib/conjugate_gradient.h"
#include "lib/matrix.h"
#include <iostream>
#include <iomanip>
using namespace std;

int n;
ColVector initx;
double alpha = 0;

Matrix getS(const ColVector &x){
    Matrix S(n, n);
    int id = 0;
    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++){
            S[i][j] = x(id);
            S[j][i] = -x(id);
            id++;
        }
    return S;
}

Matrix getQ(const Matrix &S){
    return (eye(n) + S) * inv(eye(n) - S);
}

double f(const ColVector &x){
    auto S = getS(x);
    Matrix I = eye(n);
    ColVector y = solve(I - S, initx);
    ColVector z = (I + S) * y;
    double res = pow(z(0) - alpha, 2);
    for(int i = 1; i < n; i++)
        res += pow(z(i), 2);
    return res;
}

int main(){
    cout << "Input the length of vector x:" << endl;
    cin >> n;
    initx = ColVector(n);
    cout << "Input the vector x:" << endl;
    cin >> initx;
    alpha = norm(initx);
    const int m = n * (n - 1) / 2;
    ColVector ans = CG_FR_gradfree(f, zeros(m, 1), 1e-8, 0.1, 0.4);
    auto S = getS(ans);
    cout << setprecision(12);
    cout << "Matrix S:" << endl;
    cout << S << endl;
    cout << "--------------------------" << endl;
    cout << "Check Q*x:" << endl;
    cout << getQ(S) * initx << endl;
    cout << "--------------------------" << endl;

    auto z = initx;
    z(0) += sgn(initx(0)) * alpha;
    ColVector v = z / norm(z);
    auto w = v.getSubmatrix(1, n - 1, 0, 0) / v(0);
    cout << "Exact w:" << endl;
    cout << w << endl;
    cout << "where exact S = [0, w'; -w, 0]" << endl;
    return 0;
}