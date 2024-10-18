#include "lib/linprog.h"
#include <iostream>
using namespace std;

int main(){
    int n, m;
    cout << "Input Unknown numbers" << endl;
    cin >> m;
    cout << "Input c (Maximize c' * x)" << endl;
    Matrix c(m,1);
    cin >> c;
    cout << "Input Inequations numbers" << endl;
    cin >> n;
    cout << "Input A1" << endl;
    Matrix A1(n,m);
    cin >> A1;
    cout << "Input b1 (such that A1 * x <= b1)" << endl;
    Matrix b1(n,1);
    cin >> b1;
    cout << "Input Equations numbers" << endl;
    cin >> n;
    cout << "Input A2" << endl;
    Matrix A2(n,m);
    cin >> A2;
    cout << "Input b2 (such that A2 * x == b2)" << endl;
    Matrix b2(n,1);
    cin >> b2;
    Matrix sol = linprog(c,A1,b1,A2,b2);
    if(!sol.empty())
        cout << "Solution = " << sol.T() << endl << "maxval = " << value(c.T()*sol) << endl;
    return 0;
}