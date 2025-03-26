#pragma once

#include "sparse_matrix.h"
#include "matrix.h"
#include "preconditioner.h"
#include <vector>
#include <cstring>

class amgSolver{
private:
    std::vector<SparseMatrix> Ah, Ph, Rh;
    // Ah:算子  Ph:插值  Rh:限制

    void generateGrid(const int &d);
    SparseMatrix getInterpolator(const SparseMatrix &A);
    std::vector<int> getCorsetPoints(const SparseMatrix &A);

public:
    ColVector VC(const int &d, ColVector x, const ColVector &rhs) const;
    ColVector FMG(const int &d, const ColVector &rhs) const;

    void generateGrid(const SparseMatrix &A);
    ColVector solve(const ColVector &rhs, const std::string &method, const int &maxIter, const double &eps) const;
};


class AMGPreconditioner : public Preconditioner{
private:
    amgSolver solver;

public:
    AMGPreconditioner(const SparseMatrix &A){
        solver.generateGrid(A);
    }

    ColVector vmult(const ColVector &b) const{
        return solver.FMG(0, b);
    }
};

#include "avl.h"
#include <vector>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <iomanip>

using std::string;
using std::sort;
using std::vector;
using std::max;
using std::cout;
using std::endl;

const int AMG_DIRECT_N = 16;
const int AMG_SMOOTH_ITER = 3;
const double AMG_STRONG_THERESHOLD = 0.5;

vector<int> amgSolver::getCorsetPoints(const SparseMatrix &A){
    vector<int> corsetPoints;
    const int n = A.nRows();
    // influence[i]: 强影响i的点集
    // dependence[i]: 强依赖i的点集
    vector<int> dependence[n], influence[n];
    for(int i = 0; i < n; i++){
        vector<int> neighbors = A.nonzeroIndexInRow(i);
        vector<double> vals = A.nonzeroValueInRow(i);
        // 计算非对角元的最大绝对值
        double max_nondiag = 0;
        for(int idj = 0; idj < neighbors.size(); idj++){
            int j = neighbors[idj];
            if(i == j) continue;
            max_nondiag = max( max_nondiag, fabs(vals[idj]) );
        }
        for(int idj = 0; idj < neighbors.size(); idj++){
            int j = neighbors[idj];
            if(i == j || fabs(vals[idj]) < AMG_STRONG_THERESHOLD*max_nondiag) continue;
            influence[i].push_back(j);
            dependence[j].push_back(i);
        }
    }
    vector<int> lambda(n);
    for(int i = 0; i < n; i++)
        lambda[i] = dependence[i].size();
    avlTree avl;
    avl.build(lambda);

    // 染色算法
    int assigned = 0;
    bool *isAssigned = new bool[n];
    memset(isAssigned, false, sizeof(bool)*n);
    while(assigned < n){
        int i = avl.getmaxID();
        corsetPoints.push_back(i);
        isAssigned[i] = true;
        avl.setzero(i);
        assigned++;
        for(int j : dependence[i]){
            if(isAssigned[j]) continue;
            isAssigned[j] = true;
            avl.setzero(j);
            assigned++;
        }
        for(int j : dependence[i]){
            for(int k : influence[j]){
                if(isAssigned[k]) continue;
                avl.increse(k, 1);
            }
        }
    }

    delete [] isAssigned;
    return corsetPoints;
}

SparseMatrix amgSolver::getInterpolator(const SparseMatrix &A){
    vector<Triple> elements;
    vector<int> corsetPoints = getCorsetPoints(A);
    sort(corsetPoints.begin(), corsetPoints.end());
    //for(int x : corsetPoints) cout << x << ","; cout << endl;
    const int n = A.nRows();

    // rankCorsetPoint[i]: 细网格中编号为i的点在粗网格中的编号，若i不是粗网格点，用-1标记
    int *rankCorsetPoint = new int[n];
    int rankcnt = 0;
    memset(rankCorsetPoint, -1, sizeof(int)*n);
    for(const int & x : corsetPoints){
        rankCorsetPoint[x] = rankcnt++;
    }
    for(int i = 0; i < n; i++){
        vector<int> neighbors = A.nonzeroIndexInRow(i);
        vector<double> vals = A.nonzeroValueInRow(i);

        // 计算非对角元的最大绝对值
        double max_nondiag = 0, diagCoef;
        for(int idj = 0; idj < neighbors.size(); idj++){
            int j = neighbors[idj];
            if(i == j){
                diagCoef = vals[idj];
                continue;
            }
            max_nondiag = max( max_nondiag, fabs(vals[idj]) );
        }

        // 将相邻点分类：C-粗网格中的强依赖点，Ds-不在粗网格中的强依赖点，Dw-弱依赖点（不显式存储Dw）
        vector<int> C, Ds;
        for(int idj = 0; idj < neighbors.size(); idj++){
            int j = neighbors[idj];
            if(i==j) continue;
            if( fabs(vals[idj]) < AMG_STRONG_THERESHOLD * max_nondiag ){
                diagCoef += vals[idj];
            } else if(rankCorsetPoint[j] != -1){
                C.push_back(j);
            } else {
                Ds.push_back(j);
            }
        }

        if(rankCorsetPoint[i] != -1){
            // 在粗网格中
            elements.push_back(Triple(i, rankCorsetPoint[i], 1));
        } else if(Ds.empty()){
            // 不在粗网格中，但是Ds为空，防止下面的new出现size为0的问题
            for(int j : C){
                elements.push_back(Triple(i, rankCorsetPoint[j], -A(i,j)/diagCoef));
            }
        } else {
            // 不在粗网格中，根据(8.12)计算插值系数
            double *sum_amk = new double [Ds.size()];
            for(int idm = 0; idm < Ds.size(); idm++){
                int m = Ds[idm];
                sum_amk[idm] = 0;
                for(int k : C){
                    sum_amk[idm] += A(m, k);
                }
            }
            for(int j : C){
                double w = A(i,j);
                for(int idm = 0; idm < Ds.size(); idm++){
                    int m = Ds[idm];
                    if(A(m,j)){
                        w += A(i,m) * A(m,j) / sum_amk[idm];
                    }
                }
                elements.push_back(Triple(i, rankCorsetPoint[j], -w/diagCoef));
            }
        }
    }

    delete [] rankCorsetPoint;
    return SparseMatrix(n, rankcnt, elements);
}

void amgSolver::generateGrid(const int &d){
    cout << "level " << d+1 << ":  size=" << Ah[d].nCols() << "  density=" << Ah[d].density()*100 << "%" << endl;
    if(Ah[d].nCols() <= AMG_DIRECT_N){
        return;
    }
    Ph.push_back(getInterpolator(Ah[d]));
    Rh.push_back(Ph[d].T());
    Ah.push_back( Rh[d] * Ah[d] * Ph[d] );  //Garlerkin's condition
    generateGrid(d+1);
}

void amgSolver::generateGrid(const SparseMatrix &A){
    cout << "Generating AMG hierachies..." << endl;
    int timest = clock();
    Ah.push_back(A);
    generateGrid(0);
    cout << "Generated in " << std::setprecision(3) << (double)(clock()-timest)/CLOCKS_PER_SEC << "s" << endl;
    cout << std::setprecision(6);
}

ColVector amgSolver::VC(const int &d, ColVector x, const ColVector &rhs) const{
    if(d == Ah.size()-1){
        return Ah[d].LUsolve(rhs);
    }
    for(int i = 0; i < AMG_SMOOTH_ITER; i++){
        x = Ah[d].GaussSeidel(x, rhs);
    }
    x += Ph[d] * VC(d+1, zeros(Rh[d].nRows(),1), Rh[d]*(rhs-Ah[d]*x));
    for(int i = 0; i < AMG_SMOOTH_ITER; i++){
        x = Ah[d].GaussSeidel(x, rhs);
    }
    return x;
}

ColVector amgSolver::FMG(const int &d, const ColVector &rhs) const{
    if(rhs.size() <= AMG_DIRECT_N){
        return Ah[d].LUsolve(rhs);
    }
    return VC(d, Ph[d] * FMG(d+1, Rh[d]*rhs), rhs);
}

ColVector amgSolver::solve(const ColVector &rhs, const string & method, const int &maxIter, const double &eps) const{
    cout << "--------------------------------------------------------------" << endl;
    cout << "Solving..." << endl;
    int timest = clock();
    ColVector res(rhs.size()), newres;
    for(int T = 0; T < maxIter; T++){
        if(method == "V"){
            newres = VC(0, res, rhs);
        } else {
            newres = res + FMG(0, rhs-Ah[0]*res);
        }
        double relative_err = vecnorm(newres-res, 0);
        cout << "Iteration " << T+1 << ":  residual=" << vecnorm(Ah[0]*newres-rhs, 0) << "  relative_error=" << relative_err << endl;
        res = newres;
        if(relative_err < eps) break;
    }
    cout << "Solved in " << std::setprecision(3) << (double)(clock()-timest)/CLOCKS_PER_SEC << "s" << endl;
    cout << std::setprecision(6) << "--------------------------------------------------------------" << endl;
    return res;
}
