#pragma once

#include <iostream>
#include <cstring>
#include "matrix.h"

// 带状矩阵库
template<int b1, int b2 = b1>
class BandMatrix{
private:
    int n;
    double* data[b1 + b2 + 1];

public:
    BandMatrix(int _n): n(_n){
        for(int i = -b1; i <= b2; i++){
            data[i + b1] = new double[n];
            memset(data[i + b1], 0, sizeof(double) * n);
        }
    }

    double operator () (int i, int j) const{
        return data[j - i + b1][i];
    }

    double& operator () (int i, int j){
        return data[j - i + b1][i];
    }

    int nRows() const { return n; }
    int nCols() const { return n; }

    ColVector operator * (const ColVector &x) const{
        ColVector y(n);
        for(int j = -b1; j <= b2; j++)
            for(int i = std::max(0, -j); i < n - std::max(0, j); i++)
                y(i + j) += data[j + b1][i] * x(i);
        return y;
    }

    ~BandMatrix(){
        for(int i = -b1; i <= b2; i++)
            delete[] data[i + b1];
    }

    BandMatrix& operator = (const BandMatrix &rhs){
        if(this == &rhs) return *this;
        for(int i = -b1; i <= b2; i++)
            memcpy(data[i + b1], rhs.data[i + b1], sizeof(double) * n);
        return *this;
    }

    BandMatrix T() const{
        BandMatrix res(n);
        for(int i = 0; i < n; i++)
            for(int j = std::max(0, i - b1); j <= std::min(n - 1, i + b2); j++)
                res(j, i) = data[j - i + b1][i];
        return res;
    }

    friend BandMatrix tril(const BandMatrix &rhs, int k){
        BandMatrix res(rhs.n);
        for(int i = 0; i < rhs.n; i++)
            for(int j = std::max(0, i - b1); j <= std::min(rhs.n - 1, i + b2); j++)
                if(j - i <= k) res(i, j) = rhs(i, j);
        return res;
    }

    friend BandMatrix triu(const BandMatrix &rhs, int k){
        BandMatrix res(rhs.n);
        for(int i = 0; i < rhs.n; i++)
            for(int j = std::max(0, i - b1); j <= std::min(rhs.n - 1, i + b2); j++)
                if(j - i >= k) res(i, j) = rhs(i, j);
        return res;
    }

    friend ColVector diag(const BandMatrix &rhs){
        ColVector res(rhs.n);
        for(int i = 0; i < rhs.n; i++)
            res(i) = rhs(i, i);
        return res;
    }

    void setdiag(const ColVector &rhs){
        for(int i = 0; i < n; i++)
            data[b1][i] = rhs(i);
    }

    BandMatrix operator + (const BandMatrix &rhs) const{
        BandMatrix res(n);
        for(int i = -b1; i <= b2; i++)
            for(int j = 0; j < n; j++)
                res.data[i + b1][j] = data[i + b1][j] + rhs.data[i + b1][j];
        return res;
    }

    friend ColVector solveLowerTriangular(const BandMatrix &L, const ColVector &b){
        ColVector x(b.size());
        for(int i = 0; i < L.n; i++){
            double sum = 0, diagonal = 0;
            for(int j = std::max(0, i - b1); j <= std::min(L.n - 1, i + b2); j++){
                if(j != i)
                    sum += L.data[j - i + b1][i] * x(j);
                else
                    diagonal = L.data[j - i + b1][i];
            }
            x(i) = (b(i) - sum) / diagonal;
        }
        return x;
    }

    friend ColVector solveUpperTriangular(const BandMatrix &U, const ColVector &b){
        ColVector x(b.size());
        for(int i = U.n - 1; i >= 0; i--){
            double sum = 0, diagonal = 0;
            for(int j = std::max(0, i - b1); j <= std::min(U.n - 1, i + b2); j++){
                if(j != i)
                    sum += U.data[j - i + b1][i] * x(j);
                else
                    diagonal = U.data[j - i + b1][i];
            }
            x(i) = (b(i) - sum) / diagonal;
        }
        return x;
    }
};

