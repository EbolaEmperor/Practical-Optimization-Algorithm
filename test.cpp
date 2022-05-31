#include "lib/matrix.h"
#include <iostream>
using namespace std;

int main(){
    int n;
    cin >> n;
    Matrix A(n,n);
    cin >> A;
    Matrix Q(n,n);
    Matrix T(n,n);
    Q[0][0] = 1;
    T[0][0] = A[0][0];
    Matrix tmp = A*Q.getSubmatrix(0,n-1,0,0)-T[0][0]*Q.getSubmatrix(0,n-1,0,0);
    T[0][1] = T[1][0] = tmp.vecnorm(2);
    Q.setSubmatrix(0,1,(1.0/T[0][1])*tmp);
    for(int i = 1; i < n-1; i++){
        Matrix q = Q.getSubmatrix(0,n-1,i,i);
        T[i][i] = value(q.T()*A*q);
        Matrix tmp = A*q-T[i][i]*q-T[i-1][i]*Q.getSubmatrix(0,n-1,i-1,i-1);
        T[i+1][i] = T[i][i+1] = tmp.vecnorm(2);
        Q.setSubmatrix(0,i+1, (1.0/T[i+1][i])*tmp);
    }
    Matrix q = Q.getSubmatrix(0,n-1,n-1,n-1);
    T[n-1][n-1] = value(q.T()*A*q);
    cout << Q << endl;
    cout << T << endl;
    cout << Q*T*Q.T()-A << endl;
}