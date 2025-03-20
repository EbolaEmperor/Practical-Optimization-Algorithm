#pragma once

#include "matrix.h"

class Preconditioner{
public:
    virtual ColVector vmult(const ColVector &b) const = 0;
};

template <typename Mat>
class SSORPreconditioner : public Preconditioner{
private:
    Mat M1, M2;
    ColVector d;
    int bandwidth;
    double omega;

public:
    SSORPreconditioner(const Mat &A, double omega = 1.0, int bandwidth = -1)
      : bandwidth(bandwidth), 
        omega(omega), 
        M1(A.nRows(), A.nCols()), 
        M2(A.nRows(), A.nCols())
    {
        Mat L = tril(A, -1);
        d = diag(A);
        M1.setdiag((1.0 / omega) * d);
        M1 = M1 + L;
        M2 = M1.T();
    }

    ColVector vmult(const ColVector &b) const{
        ColVector x = solveLowerTriangular(M1, b, bandwidth);
        for(int i = 0; i < d.size(); i++) x(i) *= d(i);
        x = solveUpperTriangular(M2, x, bandwidth);
        return x * ((2 - omega) / omega);
    }
};