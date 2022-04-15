/***************************************************************
 *
 * 这是一个矩阵运算库，为了方便以后设计算法更加简洁，特编写以用
 * 版本号：v1.0.2
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#include <iostream>
#include <cstring>
#include <cmath>

#ifndef _MATRIX_H_
#define _MATRIX_H_

class Matrix{
private:
    double *a;
public:
    int n, m;
    Matrix(){
        n = m = 0;
        a = nullptr;
    }
    Matrix(const int &_n, const int &_m){
        n = _n;
        m = _m;
        a = new double[n*m];
        memset(a, 0, sizeof(double)*(n*m));
    }
    Matrix(const Matrix &A){
        n = A.n;
        m = A.m;
        a = new double[n*m];
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                a[i*m+j] = A[i][j];
    }
    ~Matrix(){
        delete a;
    }

    Matrix & operator = (const Matrix & rhs){
        Matrix copy(rhs);
        std::swap(*this, copy);
        return *this;
    }
    Matrix & operator = (Matrix && rhs){
        std::swap(a, rhs.a);
        n = rhs.n;
        m = rhs.m;
        return *this;
    }

    const double* operator [] (const int &r) const{
        return a + r*m;
    }
    double* operator [] (const int &r){
        return a + r*m;
    }

    Matrix operator + (const Matrix &B) {
        if(n!=B.n || m!=B.m){
            std::cerr << "Matrix Error! Undefined Addition!" << std::endl;
            exit(-1);
        }
        Matrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C[i][j] = a[i*m+j] + B[i][j];
        return C;
    }

    Matrix operator - () {
        Matrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C[i][j] = -a[i*m+j];
        return C;
    }

    Matrix operator - (const Matrix &B) {
        if(n!=B.n || m!=B.m){
            std::cerr << "Matrix Error! Undefined Subtraction!" << std::endl;
            exit(-1);
        }
        Matrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C[i][j] = a[i*m+j] - B[i][j];
        return C;
    }

    friend Matrix operator * (const double &k, const Matrix &A) {
        Matrix C(A.n, A.m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C[i][j] = k * A[i][j];
        return C;
    }

    Matrix operator * (const Matrix &B) const{
        if(m!=B.n){
            std::cerr << "Matrix Error! Undefined multiplication!" << std::endl;
            exit(-1);
        }
        Matrix C(n, B.m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                for(int k = 0; k < m; k++)
                    C[i][j] += a[i*m+k] * B[k][j];
        return C;
    }

    Matrix T(){
        Matrix C(m, n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                C[j][i] = a[i*m+j];
        return C;
    }

    double vecnorm(const double &p){
        double norm = 0;
        for(int i = 0; i < n*m; i++)
            norm += pow(a[i], p);
        return pow(norm, 1.0/p);
    }

    friend std::ostream& operator << (std::ostream& out, const Matrix &A){
        for(int i = 0; i < A.n; i++)
        {
            out << "[ " << A[i][0];
            for(int j = 1; j < A.m; j++)
                out << ", " << A[i][j];
            out << " ]" << std::endl;
        }
        return out;
    }

    void swaprow(const int &r1, const int &r2){
        if(r1<0 || r1>=n || r2<0 || r2>=n){
            std::cerr << "Matrix Error! Swaprow ouof range!" << std::endl;
            exit(-1);
        }
        for(int j = 0; j < m; j++)
            std::swap(a[r1*m+j], a[r2*m+j]);
    }

    void swapcol(const int &r1, const int &r2){
        if(r1<0 || r1>=m || r2<0 || r2>=m){
            std::cerr << "Matrix Error! Swapcol ouof range!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < n; i++)
            std::swap(a[i*m+r1], a[i*m+r2]);
    }
};

Matrix hilbert(const int &n){
    Matrix H(n, n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            H[i][j] = 1.0/(i+j+1);
    return H;
}

Matrix zeros(const int &n, const int &m){
    return Matrix(n, m);
}

Matrix ones(const int &n, const int &m){
    Matrix H(n, m);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            H[i][j] = 1.0;
    return H;
}

Matrix eye(const int &n){
    Matrix H(n, n);
    for(int i = 0; i < n; i++)
        H[i][i] = 1.0;
    return H;
}

double value(const Matrix &A){
    return A[0][0];
}

#endif