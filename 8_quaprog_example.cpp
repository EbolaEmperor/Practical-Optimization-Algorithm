#include "lib/matrix.h"
#include "lib/quadric_programming.h"
#include <iostream>
using namespace std;

void example1(){
    const double Hv[] = {2,0,0,2};
    const double cv[] = {-2,-4};
    const double Aev[] = {};
    const double bev[] = {};
    const double Aiv[] = {-1,-1,1,0,0,1};
    const double biv[] = {-1,0,0};
    Matrix H(2,2,Hv), Ae, Ai(3,2,Aiv);
    ColVector c(2,cv), be, bi(3,biv), x0(2);
    Matrix sol = quaprog(H,c,Ae,be,Ai,bi,x0);
}

void example2(){
    const double Hv[] = {3,-1,2,-1,2,0,2,0,4};
    const double cv[] = {1,-3,-2};
    const double Aev[] = {};
    const double bev[] = {};
    const double Aiv[] = {-3,2,-5,2,-3,-2,1,0,0,0,1,0,0,0,1};
    const double biv[] = {-4,-3,0,0,0};
    Matrix H(3,3,Hv), Ae, Ai(5,3,Aiv);
    ColVector c(3,cv), be, bi(5,biv), x0(3);
    Matrix sol = quaprog(H,c,Ae,be,Ai,bi,x0);
}

int main(){
    example1();
    example2();
    return 0;
}