#ifndef _QUADRIC_PROGRAMMING_H_
#define _QUADRIC_PROGRAMMING_H_

#include "matrix.h"
#include <iostream>

/******************************************************************
 * 求解等式约束的凸二次规划问题，使用拉格朗日乘子法。问题表述：
 * min.  1/2 x^THx + c^Tx
 * s.t.  Ax=b
 * 调用方法：sol=quaprog_equ(H,c,A,b)
 * ****************************************************************/
Matrix quaprog_equ(const Matrix &H, const Matrix &c, const Matrix &A, const Matrix &b){
    Matrix IH = H.inv();
    Matrix AHA = A*IH*A.T();
    Matrix IAHA = AHA.inv();
    Matrix AIH = A*IH;
    Matrix G = IH-AIH.T()*IAHA*AIH;
    Matrix B = IAHA*AIH;
    return B.T()*b-G*c;
}

#endif