#include "lib/trust_region.h"
#include "func/rosenbrock.h"
#include <cmath>
using namespace std;
using namespace rosenbrock;

// 以扩展Rosenbrock函数为例
//********************************************************************************

int main(){
    cin >> n;
    Matrix x(n,1), y;
    for(int i = 0; i < n; i++)
        x[i][0] = (i&1) ? 1 : -1.2;
    y = trust_region(f, grad, hessian, x, 1, 0.01, 0.2, 1e-10);
    cout << "min f = f(" << y.T() << ") = " << f(y) << endl;
    cout << "||grad|| = " << grad(y).vecnorm(2) << endl;
    return 0;
}