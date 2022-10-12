/***************************************************************
 *
 * 这是一个矩阵运算库，为了方便以后设计算法更加简洁，特编写以用
 * 版本号：v1.0.2
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cstring>
#include <cmath>

class Matrix;
class ColVector;
class RowVector;

class Matrix{
private:
    double *a;
public:
    int n, m;
    Matrix();
    Matrix(const int &_n);
    Matrix(const int &_n, const int &_m);
    Matrix(const Matrix &A);
    Matrix(const double *p, const int &_n);
    Matrix(const int &_n, const int &_m, const double *p);
    ~Matrix();
    bool empty() const;

    Matrix & operator = (const Matrix & rhs);
    Matrix & operator = (Matrix && rhs);

    const double* operator [] (const int &r) const;
    double* operator [] (const int &r);
    const double element(const int &r, const int &c) const;
    double& element(const int &r, const int &c);

    friend Matrix diag(const Matrix &A);

    void setSubmatrix(const int &r, const int &c, const Matrix &rhs);
    Matrix getSubmatrix(const int &r1, const int &r2, const int &c1, const int &c2) const;

    RowVector getRow(const int &r) const;
    ColVector getCol(const int &c) const;

    Matrix operator + (const Matrix &B) const;
    Matrix operator - () const;
    Matrix operator - (const Matrix &B) const;
    Matrix operator * (const Matrix &B) const;
    Matrix T() const;

    double vecnorm(const double &p);
    void swaprow(const int &r1, const int &r2);
    void swapcol(const int &r1, const int &r2);

    double det() const;
    Matrix inv() const;
    Matrix rref() const;
    void FGdecompose(Matrix &F, Matrix &G) const;
    Matrix pinv() const;
    double sqrsum() const;
};

class RowVector: public Matrix{
public:
    RowVector(): Matrix() {};
    RowVector(const int &n): Matrix(1,n) {};
    RowVector(const int &n, const double *p): Matrix(1,n,p) {};
    RowVector(const Matrix &rhs);
    int size() const;
    const double operator [](const int &x) const;
    double & operator [] (const int &x);
    RowVector operator + (const RowVector &rhs) const;
    RowVector operator - (const RowVector &rhs) const;
    RowVector operator - () const;
    ColVector T() const;
};

class ColVector: public Matrix{
public:
    ColVector(): Matrix() {};
    ColVector(const int &n): Matrix(n,1) {};
    ColVector(const int &n, const double *p): Matrix(n,1,p) {};
    ColVector(const Matrix &rhs);
    int size() const;
    const double operator [](const int &x) const;
    double & operator [] (const int &x);
    ColVector operator + (const ColVector &rhs) const;
    ColVector operator - (const ColVector &rhs) const;
    ColVector operator - () const;
    RowVector T() const;
};

Matrix hilbert(const int &n);
Matrix zeros(const int &n, const int &m);
Matrix ones(const int &n, const int &m);
Matrix eye(const int &n);
double value(const Matrix &A);
ColVector zeroCol(const int &n);
RowVector zeroRow(const int &n);

//----------------------Matrix相关函数---------------------------

Matrix::Matrix(){
    n = m = 0;
    a = nullptr;
}
Matrix::Matrix(const int &_n){
    n = m = _n;
    a = new double[n*m];
    memset(a, 0, sizeof(double)*(n*m));
}
Matrix::Matrix(const int &_n, const int &_m){
    n = _n;
    m = _m;
    a = new double[n*m];
    memset(a, 0, sizeof(double)*(n*m));
}
Matrix::Matrix(const Matrix &A){
    n = A.n;
    m = A.m;
    a = new double[n*m];
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            a[i*m+j] = A[i][j];
}
Matrix::Matrix(const double *p, const int &_n){
    n = _n; m = 1;
    a = new double[n];
    for(int i = 0; i < n; i++)
        a[i] = p[i];
}
Matrix::Matrix(const int &_n, const int &_m, const double *p){
    n = _n; m = _m;
    a = new double[n*m];
    for(int i = 0; i < n*m; i++)
        a[i] = p[i];
}
Matrix::~Matrix(){
    delete [] a;
}
bool Matrix::empty() const{
    return n==0 || m==0;
}

Matrix & Matrix::operator = (const Matrix & rhs){
    Matrix copy(rhs);
    std::swap(*this, copy);
    return *this;
}
Matrix & Matrix::operator = (Matrix && rhs){
    std::swap(a, rhs.a);
    n = rhs.n;
    m = rhs.m;
    return *this;
}

const double* Matrix::operator [] (const int &r) const{
    return a + r*m;
}
double* Matrix::operator [] (const int &r){
    return a + r*m;
}
const double Matrix::element(const int &r, const int &c) const{
    return a[r*m+c];
}
double& Matrix::element(const int &r, const int &c){
    return a[r*m+c];
}

Matrix diag(const Matrix &A){
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
void Matrix::setSubmatrix(const int &r, const int &c, const Matrix &rhs){
    if(r<0 || c<0 || r+rhs.n>n || c+rhs.m>m){
        std::cerr << "Matrix Error! setSubmatrix::: out of range!" << std::endl;
        exit(-1);
    }
    for(int i = 0; i < rhs.n; i++)
        for(int j = 0; j < rhs.m; j++)
            a[(r+i)*m+(j+c)] = rhs[i][j];
}

Matrix Matrix::getSubmatrix(const int &r1, const int &r2, const int &c1, const int &c2) const{
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

RowVector Matrix::getRow(const int &r) const{
    RowVector res(m);
    for(int i = 0; i < m; i++)
        res[i] = a[r*m+i];
    return res;
}

ColVector Matrix::getCol(const int &c) const{
    ColVector res(n);
    for(int i = 0; i < n; i++)
        res[i] = a[i*m+c];
    return res;
}

Matrix Matrix::operator + (const Matrix &B) const {
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

Matrix Matrix::operator - () const {
    Matrix C(n, m);
    for(int i = 0; i < C.n; i++)
        for(int j = 0; j < C.m; j++)
            C[i][j] = -a[i*m+j];
    return C;
}

Matrix Matrix::operator - (const Matrix &B) const {
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

Matrix operator * (const double &k, const Matrix &A) {
    Matrix C(A.n, A.m);
    for(int i = 0; i < C.n; i++)
        for(int j = 0; j < C.m; j++)
            C[i][j] = k * A[i][j];
    return C;
}

Matrix Matrix::operator * (const Matrix &B) const{
    if(m!=B.n){
        std::cerr << "Matrix Error! Undefined multiplication! (" << n << "*" << m << ") * (" << B.n << "*" << B.m << ")" << std::endl;
        exit(-1);
    }
    Matrix C(n, B.m);
    for(int i = 0; i < C.n; i++)
        for(int j = 0; j < C.m; j++)
            for(int k = 0; k < m; k++)
                C[i][j] += a[i*m+k] * B[k][j];
    return C;
}

Matrix Matrix::T() const{
    Matrix C(m, n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            C[j][i] = a[i*m+j];
    return C;
}

double Matrix::vecnorm(const double &p){
    double norm = 0;
    for(int i = 0; i < n*m; i++)
        norm += pow(a[i], p);
    return pow(norm, 1.0/p);
}

// 将矩阵中所有元素取绝对值后返回，用法：B=abs(A)
Matrix abs(const Matrix &A){
    Matrix C = A;
    for(int i = 0; i < A.n; i++)
        for(int j = 0; j < A.m; j++)
            C[i][j] = std::fabs(A[i][j]);
    return C;
}

// 矩阵求最大值函数，返回矩阵中的最大元素，用法：x=max(A)
double max(const Matrix &A){
    if(A.m==0 || A.n==0) return 0;
    double res = A[0][0];
    for(int i = 0; i < A.n; i++)
        for(int j = 0; j < A.m; j++)
            res = std::max(res, A[i][j]);
    return res;
}

// 矩阵求和函数，返回所有元素的和，用法：s=sum(A)
double sum(const Matrix &A){
    double res = 0;
    for(int i = 0; i < A.n; i++)
        for(int j = 0; j < A.m; j++)
            res += A[i][j];
    return res;
}

std::istream& operator >> (std::istream& in, Matrix &A){
    for(int i = 0; i < A.n; i++)
        for(int j = 0; j < A.m; j++)
            in >> A[i][j];
    return in;
}

std::ostream& operator << (std::ostream& out, const Matrix &A){
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

void Matrix::swaprow(const int &r1, const int &r2){
    if(r1<0 || r1>=n || r2<0 || r2>=n){
        std::cerr << "Matrix Error! Swaprow ouof range!" << std::endl;
        exit(-1);
    }
    for(int j = 0; j < m; j++)
        std::swap(a[r1*m+j], a[r2*m+j]);
}

void Matrix::swapcol(const int &r1, const int &r2){
    if(r1<0 || r1>=m || r2<0 || r2>=m){
        std::cerr << "Matrix Error! Swapcol ouof range!" << std::endl;
        exit(-1);
    }
    for(int i = 0; i < n; i++)
        std::swap(a[i*m+r1], a[i*m+r2]);
}

// 向量对应元素相除除，用法：x = dotdiv(a,b)
Matrix dotdiv(const Matrix &a, const Matrix &b){
    Matrix x = a;
    for(int i = 0; i < a.n; i++)
        for(int j = 0; j < a.m; j++)
            x[i][j] /= b[i][j];
    return x;
}

// 解下三角方程，用法：x = solveLowerTriangular(A,b)
Matrix solveLowerTriangular(const Matrix &A, const Matrix &b){
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
Matrix solveUpperTriangular(const Matrix &A, const Matrix &b){
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
Matrix solve(Matrix A, Matrix b){
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
        if(p!=i) A.swaprow(i,p), b.swaprow(i,p);
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

// 返回矩阵的行列式，用高斯消元法计算。用法：d=A.det()
double Matrix::det() const{
    if(m!=n || m==0){
        std::cerr << "fracMatrix Error! Cannot calculate the determinate of a non-square or empty matrix!" << std::endl;
        return 0;
    }
    Matrix A(*this);
    double ans = 1;
    for(int i = 0; i < n; i++){
        int p = i;
        while(p<n && A[p][i]==0) p++;
        if(p==n) return 0;
        if(p!=i) A.swaprow(i,p);
        ans *= A[i][i];
        for(int j = i+1; j < n; j++){
            double coef = A[j][i]/A[i][i];
            for(int k = i; k < n; k++)
                A[j][k] -= A[i][k]*coef;
        }
    }
    return ans;
}

// 返回矩阵的行列式，只是提供A.det()方法的另一种调用方式。用法：d=det(A)
double det(const Matrix &A){
    return A.det();
}

// 返回矩阵的逆矩阵，用列主元高斯消元法计算。用法：B=A.inv()
Matrix Matrix::inv() const{
    if(m!=n || m==0){
        std::cerr << "Matrix Error! There's no inverse of a non-square or empty matrix!" << std::endl;
        return Matrix();
    }
    Matrix A(n,2*n);
    for(int i = 0; i < n; i++)
    {
        A[i][i+n] = 1;
        for(int j = 0; j < n; j++)
            A[i][j] = a[i*m+j];
    }
    for(int i = 0; i < n; i++){
        int p = i;
        for(int j = i+1; j < n; j++)
            if(fabs(A[j][i])>fabs(A[p][i])) p=j;
        if(A[p][i]==0){
            std::cerr << "Matrix Error! There's no inverse of a singular matrix!" << std::endl;
            return Matrix();
        }
        if(p!=i) A.swaprow(i,p);
        for(int k = i+1; k < 2*n; k++)
            A[i][k] /= A[i][i];
        A[i][i] = 1;
        for(int j = 0; j < n; j++){
            if(i==j) continue;
            double coef = A[j][i];
            for(int k = i; k < 2*n; k++)
                A[j][k] -= A[i][k]*coef;
        }
    }
    Matrix ans(n,n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            ans[i][j] = A[i][j+n];
    return ans;
}

// 返回矩阵的逆矩阵，只是提供A.inv()方法的另一种调用方式。用法：B=inv(A)
Matrix inv(const Matrix &A){
    return A.inv();
}

// 改进Cholesky分解（LDL分解），用法：L=choleskyImproved(A)，D的元素存储在L的对角线上
Matrix choleskyImproved(const Matrix &A){
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
Matrix solveByLDL(const Matrix &A, const Matrix &b){
    Matrix L = choleskyImproved(A);
    Matrix D = diag(L);
    L = L - diag(D) + eye(L.n);
    Matrix y = solveLowerTriangular(L, b);
    y = dotdiv(y,D);
    return solveUpperTriangular(L.T(), y);
}

// Gill-Murray修正Cholesky分解
Matrix gillMurray(Matrix A){
    if(A.m!=A.n || A.n==0){
        std::cerr << "Matrix Error! The method gillMurray() cannot apply on a non-square or empty matrix!" << std::endl;
        exit(-1);
    }
    using std::max;
    int n = A.n;
    double gamma=0, xi=0;
    for(int i = 0; i < n; i++)
    {
        gamma = max(gamma, fabs(A[i][i]));
        for(int j = 0; j < i; j++)
            xi = max(xi, fabs(A[i][j]));
    }
    double nu = max(1.0, sqrt(n*n-1));
    double beta2 = max(max(gamma,xi/nu),1e-6);
    Matrix c(n,n);
    for(int i = 0; i < n; i++)
        c[i][i] = A[i][i];
    Matrix L(n,n);
    for(int j = 0; j < n; j++)
    {
        int q = j;
        for(int k = j+1; k < n; k++)
            if(fabs(c[k][k])>fabs(c[q][q])) q = k;
        if(j!=q) A.swapcol(j,q), A.swaprow(j,q);
        for(int k = 0; k < j; k++)
            L[j][k] = c[j][k]/L[k][k];
        for(int i = j+1; i < n; i++)
        {
            c[i][j] = A[i][j];
            for(int k = 0; k < j; k++)
                c[i][j] -= c[i][k]*L[j][k];
        }
        double theta = 0;
        for(int k = j+1; k < n; k++)
            theta = max(theta, fabs(c[k][j]));
        L[j][j] = max(max(fabs(c[j][j]),theta*theta/beta2),1e-3);
        for(int i = j+1; i < n; i++)
            c[i][i] -= c[i][j]*c[i][j]/L[j][j];
    }
    return L;
}

// 用Gill-Murray修正Cholesky方法求解正定对称方程，用法：x = solveByLDL_GM(A,b)
Matrix solveByLDL_GM(const Matrix &A, const Matrix &b){
    Matrix L = gillMurray(A);
    Matrix D = diag(L);
    L = L - diag(D) + eye(L.n);
    Matrix y = solveLowerTriangular(L, b);
    y = dotdiv(y,D);
    return solveUpperTriangular(L.T(), y);
}

// 返回矩阵的行最简形式，调用方式：H=A.rref()
Matrix Matrix::rref() const{
    Matrix A(*this);
    int row = 0;
    for(int i = 0; i < std::min(n,m); i++){
        int p = row;
        for(int j = row+1; j < n; j++)
            if(fabs(A[j][i])>fabs(A[p][i])) p=j;
        if(A[p][i]==0) continue;
        if(p!=row) A.swaprow(row,p);
        for(int k = i+1; k < m; k++)
            A[row][k] /= A[row][i];
        A[row][i] = 1;
        for(int j = 0; j < n; j++){
            if(j==row) continue;
            double coef = A[j][i]/A[row][i];
            for(int k = i; k < m; k++)
                A[j][k] -= A[row][k]*coef;
        }
        row++;
    }
    int r = 0;
    while(r<n && A.getSubmatrix(r,r,0,m-1).vecnorm(2)>1e-10) r++;
    return r ? A.getSubmatrix(0,r-1, 0, m-1) : Matrix();
}

// 求解满秩分解A=FG，用法：A.FGdecompose(F,G)
void Matrix::FGdecompose(Matrix &F, Matrix &G) const{
    G = rref();
    F = Matrix(n, G.n);
    int num = 0;
    for(int i = 0; i < G.n; i++)
        for(int j = i; j < G.m; j++)
            if(G[i][j]==1){
                F.setSubmatrix(0,num,getSubmatrix(0,n-1,j,j));
                num++;
                break;
            }
}

// 求解矩阵的M-P广义逆，调用方法：B=A.pinv()
Matrix Matrix::pinv() const{
    Matrix F,G;
    FGdecompose(F,G);
    return G.T()*(G*G.T()).inv()*(F.T()*F).inv()*F.T();
}

// 求解矩阵M-P的广义逆，调用方法：B=pinv(A)，与B=A.pinv()等价
Matrix pinv(const Matrix &A){
    return A.pinv();
}

Matrix mergeCol(const Matrix &A, const Matrix &B){
    if(A.m == 0) return B;
    if(B.m == 0) return A;
    if(A.n != B.n){
        std::cerr << "Cannot merge Matrices in column with different row size!" << std::endl;
        return Matrix();
    }
    Matrix C(A.n, A.m+B.m);
    C.setSubmatrix(0,0,A);
    C.setSubmatrix(0,A.m,B);
    return C;
}

Matrix mergeRow(const Matrix &A, const Matrix &B){
    if(A.n == 0) return B;
    if(B.n == 0) return A;
    if(A.m != B.m){
        std::cerr << "Cannot merge Matrices in row with different column size!" << std::endl;
        return Matrix();
    }
    Matrix C(A.n+B.n, A.m);
    C.setSubmatrix(0,0,A);
    C.setSubmatrix(A.n,0,B);
    return C;
}

// 返回矩阵中所有元素的平方和
double Matrix::sqrsum() const{
    double res = 0;
    for(int i = 0; i < n*m; i++)
        res += a[i]*a[i];
    return res;
}

// 返回一个矩阵C，使得C(i,j)=min(A(i,j),B(i,j))
Matrix min(const Matrix &A, const Matrix &B){
    if(A.n!=B.n || A.m!=B.m){
        std::cerr << "Cannot use min in two matrices of different size." << std::endl;
        return Matrix();
    }
    Matrix C(A.n,A.m);
    for(int i = 0; i < C.n; i++)
        for(int j = 0; j < C.m; j++)
            C[i][j] = std::min(A[i][j],B[i][j]);
    return C;
}

// 返回一个矩阵C，使得C(i,j)=max(A(i,j),B(i,j))
Matrix max(const Matrix &A, const Matrix &B){
    if(A.n!=B.n || A.m!=B.m){
        std::cerr << "Cannot use max in two matrices of different size." << std::endl;
        return Matrix();
    }
    Matrix C(A.n,A.m);
    for(int i = 0; i < C.n; i++)
        for(int j = 0; j < C.m; j++)
            C[i][j] = std::max(A[i][j],B[i][j]);
    return C;
}

//----------------------ColVector相关函数---------------------------

ColVector::ColVector(const Matrix &rhs){
    if(rhs.m != 1){
        std::cerr << "Error! Cannot convert a non sigle-column matrix into ColVector!" << std::endl;
        (*this) = ColVector();
    } else {
        (*this) = ColVector(rhs.n);
        for(int i = 0; i < n; i++)
            (*this)[i] = rhs[i][0];
    }
}

const double ColVector::operator [](const int &x) const{
    return element(x,0);
}
    
double & ColVector::operator [] (const int &x){
    return element(x,0);
}

int ColVector::size() const{
    return n;
}

ColVector ColVector::operator + (const ColVector &rhs) const{
    if(n!=rhs.n){
        std::cerr << "ColVector Addition Error!" << std::endl;
        return ColVector();
    }
    ColVector res(n);
    for(int i = 0; i < n; i++)
        res[i] = (*this)[i] + rhs[i];
    return res;
}

ColVector ColVector::operator - (const ColVector &rhs) const{
    if(n!=rhs.n){
        std::cerr << "ColVector Addition Error!" << std::endl;
        return ColVector();
    }
    ColVector res(n);
    for(int i = 0; i < n; i++)
        res[i] = (*this)[i] - rhs[i];
    return res;
}

ColVector ColVector::operator - () const{
    ColVector res(n);
    for(int i = 0; i < n; i++)
        res[i] = -(*this)[i];
    return res;
}

ColVector operator * (const double &k, const ColVector &x){
    ColVector res(x.n);
    for(int i = 0; i < res.n; i++)
        res[i] = k * x[i];
    return res;
}

ColVector operator * (const Matrix &A, const ColVector &x){
    if(A.m!=x.n){
        std::cerr << "ColVector Multiplication Error!" << std::endl;
        return ColVector();
    }
    ColVector res(A.n);
    for(int i = 0; i < res.n; i++)
        for(int j = 0; j < A.m; j++)
            res[i] += A[i][j] * x[j];
    return res;
}

RowVector ColVector::T() const{
    RowVector res(n);
    for(int i = 0; i < n; i++)
        res[i] = (*this)[i];
    return res;
}


//----------------------RowVector相关函数---------------------------

RowVector::RowVector(const Matrix &rhs){
    if(rhs.n != 1){
        std::cerr << "Error! Cannot convert a non sigle-row matrix into RowVector!" << std::endl;
        (*this) = RowVector();
    } else {
        (*this) = RowVector(rhs.m);
        for(int i = 0; i < m; i++)
            (*this)[i] = rhs[0][i];
    }
}

const double RowVector::operator [](const int &x) const{
    return element(0,x);
}

double & RowVector::operator [] (const int &x){
    return element(0,x);
}

int RowVector::size() const{
    return m;
}

RowVector operator * (const double &k, const RowVector &x){
    RowVector res(x.n);
    for(int i = 0; i < res.n; i++)
        res[i] = k * x[i];
    return res;
}

RowVector operator * (const RowVector &x, const Matrix &A){
    if(A.n!=x.m){
        std::cerr << "RowVector Multiplication Error!" << std::endl;
        return RowVector();
    }
    RowVector res(A.m);
    for(int i = 0; i < res.m; i++)
        for(int j = 0; j < A.n; j++)
            res[i] += A[j][i] * x[j];
    return res;
}

double operator * (const RowVector &r, const ColVector &c){
    if(r.m!=c.n){
        std::cerr << "RowVector and ColVector Multiplication Error!" << std::endl;
        return 0;
    }
    return value((Matrix)r*c);
}

RowVector RowVector::operator + (const RowVector &rhs) const{
    if(m!=rhs.m){
        std::cerr << "RowVector Addition Error!" << std::endl;
        return RowVector();
    }
    RowVector res(m);
    for(int i = 0; i < m; i++)
        res[i] = (*this)[i] + rhs[i];
    return res;
}

RowVector RowVector::operator - (const RowVector &rhs) const{
    if(m!=rhs.m){
        std::cerr << "RowVector Addition Error!" << std::endl;
        return RowVector();
    }
    RowVector res(m);
    for(int i = 0; i < m; i++)
        res[i] = (*this)[i] - rhs[i];
    return res;
}

RowVector RowVector::operator - () const{
    RowVector res(m);
    for(int i = 0; i < m; i++)
        res[i] = -(*this)[i];
    return res;
}

ColVector RowVector::T() const{
    ColVector res(n);
    for(int i = 0; i < n; i++)
        res[i] = (*this)[i];
    return res;
}

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

ColVector zeroCol(const int &n){
    return ColVector(n);
}

RowVector zeroRow(const int &n){
    return RowVector(n);
}

int sgn(const double &x){
    if(x>0) return 1;
    else if(x<0) return -1;
    else return 0;
}

#endif