#include <iostream>
#include "lib/conjugate_gradient_accurate.h"
using namespace std;

int main(){
    int n;
    cin >> n;
    fracMatrix x = CG_accurate(hilbert(n), ones(n,1), zeros(n,1));
    cout << "ans = " << x.T() << endl;
    return 0;
}