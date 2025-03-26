#pragma once

#include "matrix.h"
#include <vector>

struct SparseElement{
    int j;
    double value;
    SparseElement(const int &j=0, const double &value=0):
        j(j), value(value){}
};

struct Triple{
    int i, j;
    double value;
    Triple(): i(0), j(0), value(0.0) {}
    Triple(const int &i, const int &j, const double &value):
        i(i), j(j), value(value){}
    ~Triple(){}
    bool operator < (const Triple &b) const{
        return i==b.i && j<b.j || i<b.i;
    }
};

std::vector<Triple> makeDiagTripleVector(const int &n){
    std::vector<Triple> vec;
    for(int i = 0; i < n; i++)
        vec.push_back(Triple(i, i, 0));
    return vec;
}

class SparseMatrix{
private:
    int n, m, size;
    int * row_index;
    SparseElement *elements;

public:
    int nRows() const { return n; }
    int nCols() const { return m; }

public:
    SparseMatrix();
    SparseMatrix(const int &_n, const int &_m, std::vector<Triple> ele);
    SparseMatrix(const int &_n, const int &_m)
        : SparseMatrix(_n, _m, makeDiagTripleVector(_n)) {}
    SparseMatrix(const int &_n)
        : SparseMatrix(_n, _n, makeDiagTripleVector(_n)) {}
    SparseMatrix(const SparseMatrix & rhs);
    ~SparseMatrix();

    SparseMatrix & operator = (SparseMatrix && rhs);
    SparseMatrix & operator = (const SparseMatrix & rhs);

    void clear();
    ColVector operator * (const ColVector &rhs) const;
    friend RowVector operator * (const RowVector &lhs, const SparseMatrix &A);
    SparseMatrix operator + (const SparseMatrix &rhs) const;
    SparseMatrix operator - (const SparseMatrix &rhs) const;
    ColVector wJacobi(const ColVector & x, const ColVector & b, const double &w) const;
    friend std::ostream & operator << (std::ostream & out, const SparseMatrix &A);
    Matrix toDense() const;
    ColVector LUsolve(const ColVector &b) const;

    void setdiag(const ColVector &d);
    friend ColVector diag(const SparseMatrix &A);
    friend SparseMatrix tril(const SparseMatrix &A, const int &k);

    SparseMatrix T() const;

    friend ColVector solveLowerTriangular(const SparseMatrix &L, const ColVector &b, int bandwidth);
    friend ColVector solveUpperTriangular(const SparseMatrix &U, const ColVector &b, int bandwidth);

    ColVector GaussSeidel(ColVector x, const ColVector & b) const;

    SparseMatrix operator * (const SparseMatrix &rhs) const;
    std::vector<int> nonzeroIndexInRow(const int &r) const;
    std::vector<double> nonzeroValueInRow(const int &r) const;
    double density() const;

    double operator () (const int &i, const int &j) const;
};

#include <vector>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <unordered_map>
using namespace std;

SparseMatrix::SparseMatrix(){
    n = m = size = 0;
    row_index = nullptr;
    elements = nullptr;
}

SparseMatrix::SparseMatrix(const int &_n, const int &_m, vector<Triple> ele){
    n = _n;
    m = _m;
    row_index = new int[n+1];
    sort(ele.begin(), ele.end());
    int num = 0;
    for(int i = 0; i < ele.size(); i++){
        num++;
        if(i && ele[i].i==ele[i-1].i && ele[i].j==ele[i-1].j)
            num--;
    }
    elements = new SparseElement[num];
    size = num;
    num = 0;
    int row = 0;
    row_index[0] = 0;
    for(int i = 0; i < ele.size(); i++){
        if(i && ele[i].i==ele[i-1].i && ele[i].j==ele[i-1].j){
            elements[num].value += ele[i].value;
        } else {
            while(ele[i].i > row)
                row_index[++row] = num;
            elements[num] = SparseElement(ele[i].j, ele[i].value);
            num++;
        }
    }
    while(row < n)
        row_index[++row] = size;
}

SparseMatrix::SparseMatrix(const SparseMatrix & rhs){
    n = rhs.n;
    m = rhs.m;
    size = rhs.size;
    row_index = new int[n+1];
    memcpy(row_index, rhs.row_index, sizeof(int)*(n+1));
    elements = new SparseElement[size];
    memcpy(elements, rhs.elements, sizeof(SparseElement)*size);
}

SparseMatrix& SparseMatrix::operator = (SparseMatrix && rhs){
    n = rhs.n;
    m = rhs.m;
    size = rhs.size;
    row_index = rhs.row_index;
    elements = rhs.elements;
    rhs.n = rhs.m = rhs.size = 0;
    rhs.row_index = nullptr;
    rhs.elements = nullptr;
    return *this;
}

SparseMatrix& SparseMatrix::operator = (const SparseMatrix & rhs){
    if(this == &rhs) return *this;
    clear();
    n = rhs.n;
    m = rhs.m;
    size = rhs.size;
    row_index = new int[n+1];
    memcpy(row_index, rhs.row_index, sizeof(int)*(n+1));
    elements = new SparseElement[size];
    memcpy(elements, rhs.elements, sizeof(SparseElement)*size);
    return *this;
}

SparseMatrix::~SparseMatrix(){
    clear();
}

void SparseMatrix::clear(){
    n = m = size = 0;
    delete [] row_index;
    delete [] elements;
    row_index = nullptr;
    elements = nullptr;
}

SparseMatrix SparseMatrix::operator + (const SparseMatrix &rhs) const{
    if(n != rhs.n || m != rhs.m){
        cerr << "[Error] Cannot use operator + at matrixs of distinct size!" << endl;
        exit(-1);
    }
    vector<Triple> vec;
    for(int i = 0; i < n; i++)
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            vec.push_back(Triple(i, elements[j].j, elements[j].value));
    for(int i = 0; i < n; i++)
        for(int j = rhs.row_index[i]; j < rhs.row_index[i+1]; j++)
            vec.push_back(Triple(i, rhs.elements[j].j, rhs.elements[j].value));
    return SparseMatrix(n, m, vec);
}

SparseMatrix SparseMatrix::operator - (const SparseMatrix &rhs) const{
    if(n!=rhs.n || m!=rhs.m){
        cerr << "[Error] Cannot use operator + at matrixs of distinct size!" << endl;
        exit(-1);
    }
    vector<Triple> vec;
    for(int i = 0; i < n; i++)
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            vec.push_back(Triple(i, elements[j].j, elements[j].value));
    for(int i = 0; i < n; i++)
        for(int j = rhs.row_index[i]; j < rhs.row_index[i+1]; j++)
            vec.push_back(Triple(i, rhs.elements[j].j, -rhs.elements[j].value));
    return SparseMatrix(n, m, vec);
}

ColVector SparseMatrix::operator * (const ColVector & rhs) const{
    if(nCols() != rhs.size()){
        cerr << "[Error] The columns of SparseMatrix does not coincide the rows of ColVector!" << endl;
        exit(-1);
    }
    ColVector res(nRows());
    for(int i = 0; i < res.size(); i++){
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            res(i) += elements[j].value * rhs(elements[j].j);
    }
    return res;
}

RowVector operator * (const RowVector & lhs, const SparseMatrix &A){
    if(A.nRows() != lhs.size()){
        cerr << "[Error] The rows of SparseMatrix does not coincide the columns of RowVector!" << endl;
        exit(-1);
    }
    RowVector res(A.nCols());
    for(int i = 0; i < lhs.size(); i++){
        for(int j = A.row_index[i]; j < A.row_index[i+1]; j++)
            res(A.elements[j].j) += A.elements[j].value * lhs(i);
    }
    return res;
}

ColVector SparseMatrix::wJacobi(const ColVector &x, const ColVector &b, const double &w) const{
    ColVector x1 = b;
    for(int i = 0; i < n; i++){
        double coef = 0;
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            if(elements[c].j!=i) x1(i) -= elements[c].value * x(elements[c].j);
            else coef = elements[c].value;
        x1(i) /= coef;
    }
    return (1-w)*x + w*x1;
}

ostream & operator << (std::ostream & out, const SparseMatrix &A){
    out << "shape: " << A.n << " * " << A.m << endl;
    out << "non-zero elements:" << endl;
    for(int i = 0; i < A.n; i++)
        for(int j = A.row_index[i]; j < A.row_index[i+1]; j++)
            out << "(" << i << ", " << A.elements[j].j << ", " << A.elements[j].value << ")"<< std::endl;
    out << "row_index:" << endl;
    for(int i = 0; i <= A.n; i++)
        out << A.row_index[i] << ", ";
    out << endl;
    return out;
}

Matrix SparseMatrix::toDense() const{
    Matrix A(n,m);
    for(int i = 0; i < n; i++)
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            A(i, elements[c].j) = elements[c].value;
    return A;
}

ColVector SparseMatrix::LUsolve(const ColVector &b) const{
    Matrix A = toDense();
    return A.solve(b);
}

ColVector diag(const SparseMatrix &A){
    ColVector res(A.n);
    for(int i = 0; i < A.n; i++)
        for(int j = A.row_index[i]; j < A.row_index[i+1]; j++)
            if(A.elements[j].j==i) res(i) = A.elements[j].value;
    return res;
}

SparseMatrix tril(const SparseMatrix &A, const int &k){
    vector<Triple> vec;
    for(int i = 0; i < A.n; i++)
        for(int j = A.row_index[i]; j < A.row_index[i+1]; j++)
            if(A.elements[j].j <= i + k)
                vec.push_back(Triple(i, A.elements[j].j, A.elements[j].value));
    return SparseMatrix(A.n, A.m, vec);
}

void SparseMatrix::setdiag(const ColVector &d){
    for(int i = 0; i < n; i++){
        bool flag = false;
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            if(elements[j].j == i){
                elements[j].value = d(i);
                flag = true;
            }
        if(!flag){
            cerr << "[Error] Cannot find diagonal element in row " << i << "!" << endl;
            exit(-1);
        }
    }
}

SparseMatrix SparseMatrix::T() const{
    vector<Triple> vec;
    for(int i = 0; i < n; i++)
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            vec.push_back(Triple(elements[j].j, i, elements[j].value));
    return SparseMatrix(m, n, vec);
}

ColVector solveLowerTriangular(const SparseMatrix &L, const ColVector &b, int bandwidth = -1){
    if(L.n != L.m){
        cerr << "[Error] Cannot solve a non-square matrix!" << endl;
        exit(-1);
    }
    if(b.n != L.n){
        cerr << "[Error] The size of ColVector does not coincide the size of SparseMatrix!" << endl;
        exit(-1);
    }
    ColVector x(b.n);
    for(int i = 0; i < L.n; i++){
        double sum = 0, diagonal = 0;
        for(int j = L.row_index[i]; j < L.row_index[i+1]; j++){
            if(L.elements[j].j != i)
                sum += L.elements[j].value * x(L.elements[j].j);
            else
                diagonal = L.elements[j].value;
        }
        if(diagonal == 0){
            cerr << "[Error] The diagonal element is zero!" << endl;
            exit(-1);
        }
        x(i) = (b(i) - sum) / diagonal;
    }
    return x;
}

ColVector solveUpperTriangular(const SparseMatrix &U, const ColVector &b, int bandwidth = -1){
    if(U.n != U.m){
        cerr << "[Error] Cannot solve a non-square matrix!" << endl;
        exit(-1);
    }
    if(b.n != U.n){
        cerr << "[Error] The size of ColVector does not coincide the size of SparseMatrix!" << endl;
        exit(-1);
    }
    ColVector x(b.n);
    for(int i = U.n - 1; i >= 0; i--){
        double sum = 0, diagonal = 0;
        for(int j = U.row_index[i]; j < U.row_index[i+1]; j++){
            if(U.elements[j].j != i)
                sum += U.elements[j].value * x(U.elements[j].j);
            else
                diagonal = U.elements[j].value;
        }
        if(diagonal == 0){
            cerr << "[Error] The diagonal element is zero!" << endl;
            exit(-1);
        }
        x(i) = (b(i) - sum) / diagonal;
    }
    return x;
}

double SparseMatrix::density() const{
    return (double)size/((double)n*m);
}

SparseMatrix SparseMatrix::operator * (const SparseMatrix &rhs) const{
    unordered_map<long long, double> f;
    for(int i = 0; i < n; i++)
        for(int c = row_index[i]; c < row_index[i+1]; c++){
            int k = elements[c].j;
            for(int s = rhs.row_index[k]; s < rhs.row_index[k+1]; s++){
                f[ 1ll*i*rhs.m + rhs.elements[s].j ] += elements[c].value * rhs.elements[s].value;
            }
        }
    vector<Triple> elem;
    for(auto p : f){
    	if(fabs(p.second)<1e-16) continue;
        int r = p.first / rhs.m;
        int c = p.first - 1ll*rhs.m*r;
        elem.push_back(Triple(r,c,p.second));
    }
    return SparseMatrix(n, rhs.m, elem);
}

double SparseMatrix::operator () (const int &i, const int &j) const{
    for(int c = row_index[i]; c < row_index[i+1]; c++){
        if(elements[c].j == j){
            return elements[c].value;
        } else if(elements[c].j > j){
            break;
        }
    }
    return 0;
}

vector<int> SparseMatrix::nonzeroIndexInRow(const int &r) const{
    vector<int> p;
    for(int c = row_index[r]; c < row_index[r+1]; c++)
        p.push_back(elements[c].j);
    return p;
}

vector<double> SparseMatrix::nonzeroValueInRow(const int &r) const{
    vector<double> p;
    for(int c = row_index[r]; c < row_index[r+1]; c++)
        p.push_back(elements[c].value);
    return p;
}

ColVector SparseMatrix::GaussSeidel(ColVector x, const ColVector &b) const{
    for(int i = 0; i < n; i++){
        double coef = 0, sum = 0;
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            if(elements[c].j!=i) sum += elements[c].value * x(elements[c].j);
            else coef = elements[c].value;
        x(i) = (b(i) - sum) / coef;
    }
    return x;
}