#include "lib/matrix_fraction.h"
using namespace std;

int main() {
    vector<fraction> val = {4, -1, 0, 0,
                            4, 4, -1, 0,
                            0, 10, 4, -1,
                            0, 0, 18, 4};
    fracMatrix A(4, 4, val);
    auto [P, L, U] = LU(A);
    cout << "P = \n" << P << endl << endl;
    cout << "L = \n" << L << endl << endl;
    cout << "U = \n" << U << endl << endl;
    cout << "L*U = \n" << L*U << endl;

    vector<fraction> b_val = {1, 0, 1, 0};
    fracMatrix b(4, 1, b_val);
    auto x = solve(A, b);
    cout << "A\\b = " << x.T() << endl;
    return 0;
}