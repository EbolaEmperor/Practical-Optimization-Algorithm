#include "lib/polynomial.h"
using namespace std;

void output(Complex x){
    if(fabs(x.imag()) < 1e-12) cout << x.real();
    else if(x.imag() < 0) cout << x.real() << " - " << -x.imag() << "i";
    else cout << x.real() << " + " << x.imag() << "i";
}

// 多项式求根、求极值程序，对某些多项式可能会卡死

int main(){
    cout << "Input the order of the polynomial:" << endl;
    int n;
    cin >> n;
    double *p = new double[n + 1];
    cout << "Input the coefficients from 0-order to n-order terms:" << endl;
    for(int i = 0; i <= n; i++)
        cin >> p[i];
    Polynomial poly(n, p);
    auto roots = poly.roots();
    cout << "The roots are:" << endl;
    for(auto x : roots)
        output(x), cout << endl;
    auto dpoly = poly.derivative();
    roots = dpoly.roots();
    cout << "The extreme or saddle points are:" << endl;
    for(auto x : roots){
        cout << "p(";
        output(x);
        cout << ") = ";
        output(poly(x));
        cout << endl;
    }
    return 0;
}