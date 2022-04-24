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

class Matrix;

Matrix hilbert(const int &n);
Matrix zeros(const int &n, const int &m);
Matrix ones(const int &n, const int &m);
Matrix eye(const int &n);
double value(const Matrix &A);

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

    friend Matrix diag(const Matrix &A){
        if(A.n==1 && A.m>0){
            Matrix D(A.m, A.m);
            for(int i = 0; i < A.m; i++)
                D[i][i] = A[0][i];
            return D;
        } else if(A.m==1 && A.n>0){
            Matrix D(A.n, A.n);
            for(int i = 0; i < A.n; i++)
                D[i][i] = A[i][0];
            return D;
        } else {
            int n = std::min(A.n,A.m);
            Matrix D(n,1);
            for(int i = 0; i < n; i++)
                D[i][0] = A[i][i];
            return D;
        }
    }

    // 将矩阵以第r行第c列为左上角的子矩阵设为rhs
    void setSubmatrix(const int &r, const int &c, const Matrix &rhs){
        if(r<0 || c<0 || r+rhs.n>n || c+rhs.m>m){
            std::cerr << "Matrix Error! setSubmatrix::: out of range!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < rhs.n; i++)
            for(int j = 0; j < rhs.m; j++)
                a[(r+i)*m+(j+c)] = rhs[i][j];
    }

    Matrix getSubmatrix(const int &r1, const int &r2, const int &c1, const int &c2) const{
        if(r1<0 || c1<0 || r2>=n || c2>=m || r2<r1 || c2<c1){
            std::cerr << "Matrix Error! getSubmatrix::: out of range!" << std::endl;
            exit(-1);
        }
        Matrix sub(r2-r1+1, c2-c1+1);
        for(int i = 0; i < r2-r1+1; i++)
            for(int j = 0; j < c2-c1+1; j++)
                sub[i][j] = a[(r1+i)*m+(c1+j)];
        return sub;
    }

    Matrix operator + (const Matrix &B) const {
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

    Matrix operator - () const {
        Matrix C(n, m);
        for(int i = 0; i < C.n; i++)
            for(int j = 0; j < C.m; j++)
                C[i][j] = -a[i*m+j];
        return C;
    }

    Matrix operator - (const Matrix &B) const {
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

    Matrix T() const{
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
            out << " ]";
            if(i < A.n-1) out << std::endl;
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

    // 向量对应元素相除除，用法：x = dotdiv(a,b)
    friend Matrix dotdiv(const Matrix &a, const Matrix &b){
        Matrix x = a;
        for(int i = 0; i < a.n; i++)
            for(int j = 0; j < a.m; j++)
                x[i][j] /= b[i][j];
        return x;
    }

    // 解下三角方程，用法：x = solveLowerTriangular(A,b)
    friend Matrix solveLowerTriangular(const Matrix &A, const Matrix &b){
        Matrix x = b;
        int n = A.n;
        for(int i = 0; i < n; i++){
            x[i][0] /= A[i][i];
            for(int j = i+1; j < n; j++)
                x[j][0] -= x[i][0] * A[j][i];
        }
        return x;
    }

    // 解上三角方程，用法：x = solveUpperTriangular(A,b)
    friend Matrix solveUpperTriangular(const Matrix &A, const Matrix &b){
        Matrix x = b;
        int n = A.n;
        for(int i = n-1; i >= 0; i--){
            x[i][0] /= A[i][i];
            for(int j = 0; j < i; j++)
                x[j][0] -= x[i][0] * A[j][i];
        }
        return x;
    }

    // 解方程Ax=b，算法为列主元法Gauss消元，用法：x = solve(A,b)
    friend Matrix solve(Matrix A, Matrix b){
        if(A.m!=A.n || A.n!=b.n || A.m==0){
            std::cerr << "Matrix Error! The method solve() cannot solve an ill-posed equation!" << std::endl;
            return Matrix();
        }
        int n = A.n;
        Matrix x(n,1);
        for(int i = 0; i < n; i++){
            int p = i;
            for(int j=i; j<n; j++)
                if(fabs(A[j][i])>fabs(A[p][i])) p=j;
            if(p!=i) A.swaprow(i,p);
            if(!A[i][i]){
                std::cerr << "Matrix Error! The method solve() cannot solve an singular equation!" << std::endl;
                return Matrix();
            }
            for(int j = 0; j < n; j++){
                if(i==j) continue;
                double coef = A[j][i]/A[i][i];
                for(int k = i; k < n; k++)
                    A[j][k] -= A[i][k]*coef;
                b[j][0] -= b[i][0]*coef;
            }
        }
        for(int i = 0; i < n; i++)
            x[i][0] = b[i][0]/A[i][i];
        return x;
    }

    // 改进Cholesky分解（LDL分解），用法：L=choleskyImproved(A)，D的元素存储在L的对角线上
    friend Matrix choleskyImproved(const Matrix &A){
        if(A.m!=A.n || A.m==0){
            std::cerr << "Matrix Error! The method cholesky() cannot apply on a non-square or empty matrix!" << std::endl;
            return Matrix();
        }
        int n = A.n;
        Matrix L(n, n);
        for(int j = 0; j < n; j++)
        {
            L[j][j] = A[j][j];
            for(int k = 0; k < j; k++)
                L[j][j] -= L[j][k]*L[k][k]*L[j][k];
            for(int i = j+1; i < n; i++)
            {
                L[i][j] = A[i][j];
                for(int k = 0; k < j; k++)
                    L[i][j] -= L[i][k]*L[k][k]*L[j][k];
                L[i][j] /= L[j][j];
            }
        }
        return L;
    }

    // 用改进Cholesky方法求解正定对称方程，用法：x = solveByLDL(A,b)
    friend Matrix solveByLDL(const Matrix &A, const Matrix &b){
        Matrix L = choleskyImproved(A);
        Matrix D = diag(L);
        L = L - diag(D) + eye(L.n);
        Matrix y = solveLowerTriangular(L, b);
        y = dotdiv(y,D);
        return solveUpperTriangular(L.T(), y);
    }

    // Gill-Murray修正Cholesky分解
    void gillMurray(Matrix &L, Matrix &D, Matrix &E) const{
        if(m!=n || n==0){
            std::cerr << "Matrix Error! The method gillMurray() cannot apply on a non-square or empty matrix!" << std::endl;
            exit(-1);
        }
        Matrix G(*this);
        double gamma = 0, xi = 0;
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(i==j) gamma = std::max(gamma, fabs(a[i*m+j]));
                else xi = std::max(xi, fabs(a[i*m+j]));
        double nu = std::max(1.0, sqrt(n*n-1));
        double beta2 = std::max(std::max(gamma,xi/nu),1e-6);
        Matrix c(n,n);
        for(int i = 0; i < n; i++)
            c[i][i] = a[i*m+i];
        L = eye(n);
        D = zeros(1,n);
        E = zeros(1,n);
        for(int j = 0; j < n; j++){
            int q = j;
            for(int k = j+1; k < n; k++)
                if(fabs(c[k][k])>fabs(c[q][q])) q = k;
            if(j!=q) G.swaprow(j,q), G.swapcol(j,q);
            for(int k = 0; k < j; k++)
                L[j][k] = c[j][k]/D[0][k];
            if(j==0) c.setSubmatrix(j+1,j, G.getSubmatrix(j+1,n-1,j,j));
            else if(j<n-1) c.setSubmatrix(j+1,j, G.getSubmatrix(j+1,n-1,j,j) - c.getSubmatrix(j+1,n-1,0,j-1) * L.getSubmatrix(j,j,0,j-1).T() );
            double theta = 0;
            for(int k = j+1; k < n; k++)
                theta = std::max(theta, fabs(c[k][j]));
            D[0][j] = std::max(1e-3, std::max(fabs(c[j][j]), theta*theta/beta2));
            E[0][j] = D[0][j] - c[j][j];
            for(int i = j+1; i < n; i++)
                c[i][i] -= c[i][j]*c[i][j]/D[0][j];
        }
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