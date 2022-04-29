#include "lib/trust_region.h"
#include "func/rosenbrock.h"
#include <cmath>
using namespace std;
using namespace rosenbrock;

// 以扩展Rosenbrock函数为例
//********************************************************************************

int main(){
    cin >> n;
    Matrix y = trust_region(f, grad, hessian, initial(), 1, 0.01, 0.2, 1e-6);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl;
    cout << "||grad|| = " << grad(y).vecnorm(2) << endl;
    return 0;
}