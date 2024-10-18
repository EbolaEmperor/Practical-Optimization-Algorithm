#include "lib/quassi_newton.h"
#include "func/rosenbrock.h"
#include <cmath>
using namespace std;
using namespace rosenbrock;

// 以扩展Rosenbrock函数为例
//********************************************************************************

int main(){
    cin >> n;
    ColVector x=initial(), y;

    // 带Wolve准则的BFGS方法
    y = bfgs_gradfree(f, x, 1e-5, 0.1, 0.4);

    // 带Goldstein准则的BFGS方法
    y = bfgs_goldstein_gradfree(f, x, 1e-5, 0.2);

    // 带Wolve准则的DFP方法
    y = dfp_gradfree(f, x, 1e-5, 0.1, 0.4);

    // 带Wolve准则的Broyden方法
    y = broyden_gradfree(f, x, 0.8, 1e-5, 0.1, 0.4);
    return 0;
}