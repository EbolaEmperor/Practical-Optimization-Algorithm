#include <iostream>
#include "lib/least_square.h"
using namespace std;

ColVector F(const ColVector &x){
    ColVector res(2);
    res[0] = x[0]*x[0]*x[0] - x[1] - 1;
    res[1] = x[0]*x[0] - x[1];
    return res;
}

Matrix J(const ColVector &x){
    Matrix res(2,2);
    res[0][0] = 3*x[0]*x[0];
    res[0][1] = -1;
    res[1][0] = 2*x[0];
    res[1][1] = -1;
    return res;
}

int main(){
    ColVector x(2);
    x[0] = 0.7; x[1] = -0.2;
    ColVector sol1 = nonlinlsq_LM(F, J, x, 1e-10, 0.3, 0.6, 2000);
    ColVector sol2 = nonlinlsq_LM_gradfree(F, x, 1e-10, 0.3, 0.6, 2000);
    return 0;
}