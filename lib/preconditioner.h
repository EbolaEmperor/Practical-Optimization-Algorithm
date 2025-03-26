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
    double omega;

public:
    SSORPreconditioner(const Mat &A, double omega = 1.0)
      : omega(omega), 
        M1(A.nRows()), 
        M2(A.nRows())
    {
        Mat L = tril(A, -1);
        d = diag(A);
        M1.setdiag((1.0 / omega) * d);
        M1 = M1 + L;
        M2 = M1.T();
    }

    ColVector vmult(const ColVector &b) const{
        ColVector x = solveLowerTriangular(M1, b);
        for(int i = 0; i < d.size(); i++) x(i) *= d(i);
        x = solveUpperTriangular(M2, x);
        return x * ((2 - omega) / omega);
    }
};


template <typename Mat>
class JORPreconditioner : public Preconditioner{
private:
    ColVector d;
    double omega;

public:
    JORPreconditioner(const Mat &A, double omega = 1.0)
      : omega(omega), 
        d(diag(A)) {}

    ColVector vmult(const ColVector &b) const{
        return omega * dotdiv(b, d);
    }
};