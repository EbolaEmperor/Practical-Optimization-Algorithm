#define SILENCE
#include "lib/matrix.h"
#include "lib/quassi_newton.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <ctime>
using namespace std;

int n;
complex<double> *a;

complex<double> f(const complex<double> &z){
    complex<double> cur = 1, ans = 0;
    for(int i = 0; i <= n; i++){
        ans += a[i]*cur;
        cur *= z;
    }
    return ans;
}

double norm2(const complex<double> &z){
    return z.real()*z.real() + z.imag()*z.imag();
}

void printcomp(complex<double> z){
    z = complex<double>(round(z.real()*1e6)*1e-6, round(z.imag()*1e6)*1e-6);
    if(!z.real() && !z.imag()) cout << 0 << endl;
    else if(!z.imag()) cout << z.real() << endl;
    else if(!z.real()) cout << z.imag() << "i" << endl;
    else if(z.imag() < 0) cout << z.real() << z.imag() << "i" << endl;
    else cout << z.real() << "+" << z.imag() << "i" << endl;
}

double F(const Matrix &x){
    complex<double> fval = f(complex<double>(x[0][0],x[1][0]));
    return norm2(fval);
}

void divide(const complex<double> &z){
    a[0] /= -z;
    for(int k = 1; k < n; k++){
        a[k] = (a[k-1]-a[k])/z;
    }
}

int main(){
    srand(time(0));
    cout << "Input degree" << endl;
    cin >> n;
    a = new complex<double> [n+1];
    cout << "Input coefficients a0,...,an" << endl;
    for(int i = 0; i <= n; i++) cin >> a[i];
    for(; n; n--){
        Matrix initial(2,1), root;
        while(1){
            initial[0][0] = (double)rand()/RAND_MAX;
            initial[1][0] = (double)rand()/RAND_MAX;
            root = bfgs_gradfree(F, initial, 1e-7, 0.3, 0.6);
            if(fabs(F(root)) < 1e-5) break;
        }
        complex<double> lambda(root[0][0], root[1][0]);
        printcomp(lambda);
        divide(lambda);
    }
    return 0;
}