#pragma once

#include "preconditioner.h"
#include "sparse_matrix.h"

class Lap2DMGSolver {
private:
    Matrix getLaplacian2D(int n) const{
        Matrix A(n*n, n*n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++){
                A(i*n + j, i*n + j) = 4;
                if(i > 0) A(i*n + j, (i-1)*n + j) = -1;
                if(i < n - 1) A(i*n + j, (i+1)*n + j) = -1;
                if(j > 0) A(i*n + j, i*n + j - 1) = -1;
                if(j < n - 1) A(i*n + j, i*n + j + 1) = -1;
            }
        return A;
    }

    ColVector restriction(const ColVector &x, int n) const{
        int m = n / 2;
        ColVector y(m * m);
        for(int i = 0; i < m; i++)
            for(int j = 0; j < m; j++)
                y(i*m + j) = ( x(2*i*n + 2*j) + x(2*i*n + 2*j + 1) 
                             + x((2*i + 1)*n + 2*j) + x((2*i + 1)*n + 2*j + 1) );
        return y;
    }

    ColVector prolongation(const ColVector &x, int n) const{
        int m = n * 2;
        ColVector y(m * m);
        for(int i = 0; i < m; i++)
            for(int j = 0; j < m; j++)
                y(i * m + j) = x(i/2 * n + j/2);
        return y;
    }

    ColVector prolongationLinear(const ColVector &x, int n) const{
        ColVector y(4 * n * n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++){
                y[2*i*2*n + 2*j] += x[i*n + j] * 4 / 9;
                y[2*i*2*n + 2*j + 1] += x[i*n + j] * 4 / 9;
                y[(2*i + 1)*2*n + 2*j] += x[i*n + j] * 4 / 9;
                y[(2*i + 1)*2*n + 2*j + 1] += x[i*n + j] * 4 / 9;
                if(i > 0){
                    y[(2*i-1)*2*n + 2*j] += x[i*n + j] * 2 / 9;
                    y[(2*i-1)*2*n + 2*j + 1] += x[i*n + j] * 2 / 9;
                }
                if(i < n - 1){
                    y[(2*i+2)*2*n + 2*j] += x[i*n + j] * 2 / 9;
                    y[(2*i+2)*2*n + 2*j + 1] += x[i*n + j] * 2 / 9;
                }
                if(j > 0){
                    y[2*i*2*n + 2*j - 1] += x[i*n + j] * 2 / 9;
                    y[(2*i+1)*2*n + 2*j - 1] += x[i*n + j] * 2 / 9;
                }
                if(j < n - 1){
                    y[2*i*2*n + 2*j + 2] += x[i*n + j] * 2 / 9;
                    y[(2*i+1)*2*n + 2*j + 2] += x[i*n + j] * 2 / 9;
                }
                if(i > 0 && j > 0) y[(2*i-1)*2*n + 2*j - 1] += x[i*n + j] / 9;
                if(i > 0 && j < n - 1) y[(2*i-1)*2*n + 2*j + 2] += x[i*n + j] / 9;
                if(i < n - 1 && j > 0) y[(2*i+2)*2*n + 2*j - 1] += x[i*n + j] / 9;
                if(i < n - 1 && j < n - 1) y[(2*i+2)*2*n + 2*j + 2] += x[i*n + j] / 9;
            }
        return y;
    }

    ColVector wJacobi(const ColVector &x, const ColVector &b, double w = 0.8) const{
        int n = sqrt(x.size());
        ColVector middle_x(n * n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++){
                double sum = 0;
                if(i > 0) sum += x[(i-1)*n + j];
                if(i < n - 1) sum += x[(i+1)*n + j];
                if(j > 0) sum += x[i*n + j - 1];
                if(j < n - 1) sum += x[i*n + j + 1];
                middle_x[i*n + j] = (1 - w) * x[i*n + j] + w * (b[i*n + j] + sum) / 4;
            }
        return middle_x;
    }

    void SSOR(ColVector &x, const ColVector &b, double w = 1.0) const{
        int n = sqrt(x.size());
        ColVector middle_x(n * n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++){
                double sum = 0;
                if(i > 0) sum += middle_x[(i-1)*n + j];
                if(i < n - 1) sum += x[(i+1)*n + j];
                if(j > 0) sum += middle_x[i*n + j - 1];
                if(j < n - 1) sum += x[i*n + j + 1];
                middle_x[i*n + j] = (1 - w) * x[i*n + j] + w * (b[i*n + j] + sum) / 4;
            }
        for(int i = n - 1; i >= 0; i--)
            for(int j = n - 1; j >= 0; j--){
                double sum = 0;
                if(i > 0) sum += middle_x[(i-1)*n + j];
                if(i < n - 1) sum += x[(i+1)*n + j];
                if(j > 0) sum += middle_x[i*n + j - 1];
                if(j < n - 1) sum += x[i*n + j + 1];
                x[i*n + j] = (1 - w) * middle_x[i*n + j] + w * (b[i*n + j] + sum) / 4;
            }
    }

public:
    ColVector vmult(const ColVector &x) const{
        int n = sqrt(x.size());
        ColVector y(n * n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++){
                y[i*n + j] = 4 * x[i*n + j];
                if(i > 0) y[i*n + j] -= x[(i-1)*n + j];
                if(i < n - 1) y[i*n + j] -= x[(i+1)*n + j];
                if(j > 0) y[i*n + j] -= x[i*n + j - 1];
                if(j < n - 1) y[i*n + j] -= x[i*n + j + 1];
            }
        return y;
    }

    ColVector VCycle(ColVector x, const ColVector &b, int n) const{
        if(n <= 4){
            Matrix A = getLaplacian2D(n);
            return solve(A, b);
        }
        for(int i = 0; i < 4; i++)
            x = wJacobi(x, b);
        ColVector e2 = VCycle(zeros(n / 2 * n / 2, 1), 
                              restriction(b - vmult(x), n), 
                              n / 2);
        x += prolongation(e2, n / 2);
        for(int i = 0; i < 4; i++)
            x = wJacobi(x, b);
        return x;
    }

    ColVector FMGCycle(const ColVector &b, int n) const{
        if(n <= 4){
            Matrix A = getLaplacian2D(n);
            return solve(A, b);
        }
        ColVector x0 = FMGCycle(restriction(b, n), n / 2);
        return VCycle(prolongation(x0, n / 2), b, n);
    }
};


class Lap2DMGPreconditioner : public Preconditioner {
private:
    Lap2DMGSolver solver;
    int n;

public:
    Lap2DMGPreconditioner(int n) : n(n) {}

    ColVector vmult(const ColVector &b) const{
        return solver.FMGCycle(b, n);
    }
};