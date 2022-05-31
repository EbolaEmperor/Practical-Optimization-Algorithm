#include "lib/matrix.h"
#include "lib/quadric_programming.h"
#include <iostream>
using namespace std;

const double Hv[] = {2,0,0,2};
const double cv[] = {-2,-4};
const double Aev[] = {};
const double bev[] = {};
const double Aiv[] = {-1,-1,1,0,0,1};
const double biv[] = {-1,0,0};

int main(){
    Matrix H(2,2,Hv), c(2,1,cv), Ae, be, Ai(3,2,Aiv), bi(3,1,biv), x0(2,1);
    Matrix sol = quaprog(H,c,Ae,be,Ai,bi,x0);
    return 0;
}