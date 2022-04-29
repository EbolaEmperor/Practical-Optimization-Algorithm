#include <iostream>
#include "lib/conjugate_gradient.h"
using namespace std;

int main(){
    int n;
    cin >> n;
    Matrix x = CG(hilbert(n), ones(n,1), zeros(n,1), 1e-6);
    cout << "ans = " << x.T() << endl;
    return 0;
}