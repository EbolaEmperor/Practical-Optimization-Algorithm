#include <iostream>
#include "lib/least_square.h"
using namespace std;

Matrix F(const Matrix &x){
    Matrix res(2,1);
    res[0][0] = x[0][0]*x[0][0]*x[0][0] - x[1][0] - 1;
    res[1][0] = x[0][0]*x[0][0] - x[1][0];
    return res;
}

Matrix J(const Matrix &x){
    Matrix res(2,2);
    res[0][0] = 3*x[0][0]*x[0][0];
    res[0][1] = -1;
    res[1][0] = 2*x[0][0];
    res[1][1] = -1;
    return res;
}

int main(){
    Matrix x(2,1);
    x[0][0] = 0.70; x[1][0] = -0.2;
    Matrix y = nonlinlsq_LM(F, J, x, 1e-10, 0.3, 0.6, 2000);
    cout << "min f = f(" << y.T() << ") = " << fval(y) << endl;
    return 0;
}