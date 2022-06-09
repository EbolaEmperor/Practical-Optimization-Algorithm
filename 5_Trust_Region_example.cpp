#include "lib/trust_region.h"
#include "func/rosenbrock.h"
#include <cmath>
using namespace std;
using namespace rosenbrock;

// 以扩展Rosenbrock函数为例
//********************************************************************************

int main(){
    cin >> n;
    ColVector sol1 = trust_region(f, grad, hessian, initial(), 1, 0.01, 0.2, 1e-8);
    ColVector sol2 = trust_region_gradfree(f, initial(), 1, 0.01, 0.2, 1e-8);
    return 0;
}