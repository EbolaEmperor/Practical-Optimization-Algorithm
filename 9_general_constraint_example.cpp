#include "lib/matrix.h"
#include "lib/general_constraint.h"
#include <iostream>
using namespace std;

double f(const ColVector &x){
    return 1000-x[0]*x[0]-2*x[1]*x[1]-x[2]*x[2]-x[0]*x[1]-x[0]*x[2];
}

ColVector gradf(const ColVector &x){
    ColVector res(3);
    res[0] = -2*x[0]-x[1]-x[2];
    res[1] = -4*x[1]-x[0];
    res[2] = -2*x[2]-x[0];
    return res;
}

ColVector h(const ColVector &x){
    ColVector res(2);
    res[0] = 8*x[0]+14*x[1]+7*x[2]-56;
    res[1] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-25;
    return res;
}

Matrix Jh(const ColVector &x){
    Matrix res(2,3);
    res[0][0]=8; res[0][1]=14; res[0][2]=7;
    res[1][0]=2*x[0]; res[1][1]=2*x[1]; res[1][2]=2*x[2];
    return res;
}

ColVector g(const ColVector &x){
    return x;
}

Matrix Jg(const ColVector &x){
    return eye(3);
}

int main(){
    ColVector x0 = 2*ones(3,1);
    ColVector sol = general_constraint_optimal_PHR(f,h,g,gradf,Jh,Jg,x0);
    return 0;
}