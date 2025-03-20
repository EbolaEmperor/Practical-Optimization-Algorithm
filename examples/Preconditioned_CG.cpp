#include "lib/conjugate_gradient.h"
#include "lib/preconditioner.h"
using namespace std;

Matrix getLaplacian1D(int n){
    Matrix A(n, n);
    for(int i = 0; i < n; i++){
        A(i, i) = 2;
        if(i > 0) A(i, i - 1) = -1;
        if(i < n - 1) A(i, i + 1) = -1;
    }
    return A;
}

int main(){
    int n;
    cin >> n;
    Matrix A = getLaplacian1D(n);
    SSORPreconditioner P(A, 1.99999999, 1);
    cout << "Start PCG Iteration..." << endl;
    ColVector b = ones(n, 1);
    ColVector x = PCG(A, b, zeros(n, 1), P, 1e-9);
    cout << "Residule: " << norm(A * x - b) << endl;
    return 0;
}