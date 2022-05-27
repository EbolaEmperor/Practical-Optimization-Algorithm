#include "lib/matrix.h"
#include "lib/quadric_programming.h"
#include <iostream>
using namespace std;

const double Hv[] = {3,-1,0,-1,2,-1,0,-1,1};
const double cv[] = {1,1,1};
const double Av[] = {1,2,1};
const double bv[] = {4};

int main(){
    Matrix H(3,3,Hv), c(3,1,cv), A(1,3,Av), b(1,1,bv);
    Matrix sol = quaprog_equ(H,c,A,b);
    cout << "sol=" << sol.T() << endl;
    cout << "val=" << value(0.5*sol.T()*H*sol+c.T()*sol) << endl;
    return 0;
}