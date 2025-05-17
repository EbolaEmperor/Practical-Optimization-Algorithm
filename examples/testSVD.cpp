#include "lib/matrix.h"
using namespace std;

int main(){
    Matrix A = randMatrix(5, 5);
    cout << "A = " << A << endl << endl;
    auto [U, S, V] = A.svd();
    cout << "U = " << U << endl << endl;
    cout << "S = " << S << endl << endl;
    cout << "V = " << V << endl << endl;
    cout << "err = " << (A - U * S * V.T()) << endl;
    return 0;
}