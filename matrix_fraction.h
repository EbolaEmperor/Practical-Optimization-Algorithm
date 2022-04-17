/***************************************************************
 *
 * 这是一个利用高精度分数运算实现的精确矩阵运算库
 * 运算效率较低，但能保证绝对精确的运算，可用于理论研究
 * 依赖于：fraction.h, bigint.h
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#ifndef _MATRIX_FRACTION_H_
#define _MATRIX_FRACTION_H_

#include <iostream>
#include <cstring>
#include <cmath>
#include "fraction.h"

class fracMatrix{
private:
    std::vector<fraction> a;
public:
    int n, m;
    fracMatrix(){
        n = m = 0;
        a.clear();
    }
    fracMatrix(const int &_n, const int &_m){
        n = _n;
        m = _m;
        a.resize(n*m);
        for(int i = 0; i < n*m; i++)
            a[i] = fraction(0,1);
    }
    fracMatrix(const fracMatrix &A){
        n = A.n;
        m = A.m;
        a.resize(n*m);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                a[i*m+j] = A.element(i,j);
    }
    ~fracMatrix(){
        a.clear();
    }

    bool iszero() const{
        for(int i = 0; i < n*m; i++)
            if(a[i]!=fraction(0,1)) return false;
        return true; 
    }
    fracMatrix & operator = (const fracMatrix & rhs){
        fracMatrix copy(rhs);
        std::swap(*this, copy);
        return *this;
    }
    fracMatrix & operator = (fracMatrix && rhs){
        std::swap(a, rhs.a);
        n = rhs.n;
        m = rhs.m;
        return *this;
    }

    fraction & element(const int &i, const int &j){
        return a[i*m+j];
    }
    const fraction & element(const int &i, const int &j) const{
        return a[i*m+j];
    }

    fracMatrix operator + (const fracMatrix &B) {
        if(n!=B.n || m!=B.m){
            std::cerr << "fracMatrix Error! Undefined Addition!" << std::endl;
            exit(-1);
        }
        fracMatrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C.element(i,j) = a[i*m+j] + B.element(i,j);
        return C;
    }

    fracMatrix operator - () {
        fracMatrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C.element(i,j) = -a[i*m+j];
        return C;
    }

    fracMatrix operator - (const fracMatrix &B) {
        if(n!=B.n || m!=B.m){
            std::cerr << "fracMatrix Error! Undefined Subtraction!" << std::endl;
            exit(-1);
        }
        fracMatrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C.element(i,j) = a[i*m+j] - B.element(i,j);
        return C;
    }

    friend fracMatrix operator * (const fraction &k, const fracMatrix &A) {
        fracMatrix C(A.n, A.m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C.element(i,j) = k * A.element(i,j);
        return C;
    }

    fracMatrix operator * (const fracMatrix &B) const{
        if(m!=B.n){
            std::cerr << "fracMatrix Error! Undefined multiplication!" << std::endl;
            exit(-1);
        }
        fracMatrix C(n, B.m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                for(int k = 0; k < m; k++)
                    C.element(i,j) += a[i*m+k] * B.element(k,j);
        return C;
    }

    fracMatrix T(){
        fracMatrix C(m, n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                C.element(j,i) = a[i*m+j];
        return C;
    }

    double vecnorm(const double &p){
        double norm = 0;
        for(int i = 0; i < n*m; i++)
            norm += pow(a[i].to_double(), p);
        return pow(norm, 1.0/p);
    }

    friend std::ostream& operator << (std::ostream& out, const fracMatrix &A){
        for(int i = 0; i < A.n; i++)
        {
            out << "[ " << A.element(i,0);
            for(int j = 1; j < A.m; j++)
                out << ", " << A.element(i,j);
            out << " ]" << std::endl;
        }
        return out;
    }

    void swaprow(const int &r1, const int &r2){
        if(r1<0 || r1>=n || r2<0 || r2>=n){
            std::cerr << "fracMatrix Error! Swaprow ouof range!" << std::endl;
            exit(-1);
        }
        for(int j = 0; j < m; j++)
            std::swap(a[r1*m+j], a[r2*m+j]);
    }

    void swapcol(const int &r1, const int &r2){
        if(r1<0 || r1>=m || r2<0 || r2>=m){
            std::cerr << "fracMatrix Error! Swapcol ouof range!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < n; i++)
            std::swap(a[i*m+r1], a[i*m+r2]);
    }
};

fracMatrix hilbert(const int &n){
    fracMatrix H(n, n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            H.element(i,j) = fraction(1, i+j+1);
    return H;
}

fracMatrix zeros(const int &n, const int &m){
    return fracMatrix(n, m);
}

fracMatrix ones(const int &n, const int &m){
    fracMatrix H(n, m);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            H.element(i,j) = 1;
    return H;
}

fracMatrix eye(const int &n){
    fracMatrix H(n, n);
    for(int i = 0; i < n; i++)
        H.element(i,i) = 1;
    return H;
}

fraction value(const fracMatrix &A){
    return A.element(0,0);
}

#endif