#include "lib/conjugate_gradient.h"
#include "lib/preconditioner.h"
#include "lib/amg.h"
#include "lib/mg_lap2d.h"
#include "lib/sparse_matrix.h"
#include "lib/pnglib.h"
using namespace std;

Matrix getLaplacian1D(int n){
    Matrix A(n, n);
    for(int i = 0; i < n; i++){
        A(i, i) = 2;
        if(i > 0) A(i, i - 1) = -1;
        if(i < n - 1) A(i, i + 1) = -1;
    }
    return A;
}

SparseMatrix getLaplacian2D(int n){
    vector<Triple> vec;
    auto idx = [&](int i, int j){ return i * n + j; };
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++){
            vec.emplace_back(idx(i, j), idx(i, j), 4);
            if(i > 0) vec.emplace_back(idx(i, j), idx(i - 1, j), -1);
            if(i < n - 1) vec.emplace_back(idx(i, j), idx(i + 1, j), -1);
            if(j > 0) vec.emplace_back(idx(i, j), idx(i, j - 1), -1);
            if(j < n - 1) vec.emplace_back(idx(i, j), idx(i, j + 1), -1);
        }
    return SparseMatrix(n * n, n * n, vec);
}

double f(double x, double y){
    if(x + y < 1) return 1;
    return -1;
}

ColVector discreteF(int n){
    double h = 1.0 / (n + 1);
    ColVector fvec(n * n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            fvec(i * n + j) = f((i + 1) * h, (j + 1) * h);
    return fvec;
}

int main(int argc, char * argv[]){
    if(argc != 3){
        cout << "Usage: " << argv[0] << " <testID> <n>" << endl;
        return 0;
    }
    int testID = stoi(argv[1]);
    int n = stoi(argv[2]);
    const double eps = 1e-6;

    if(testID == 1){
        cout << "--------- 1D Laplacian ----------" << endl;
        cout << "Grid size: " << n << endl;
        cout << "\e[1;32mPreparing\e[0m preconditioner..." << endl;
        Matrix A = getLaplacian1D(n);
        SSORPreconditioner P(A, 1.99999999, 1);

        cout << "\e[1;32mRunning\e[0m PCG Iteration..." << endl;
        double h = 1.0 / (n + 1);
        ColVector b = ones(n, 1) * h;
        ColVector x = PCG(A, b, zeros(n, 1), P, eps);
        cout << "Residule: " << norm(A * x - b) << endl;
        cout << "---------------------------------" << endl;
    } 
    else if(testID == 2){
        cout << "--------- SSOR-Preconditioned CG For 2D Laplacian ----------" << endl;
        cout << "Grid size: " << n << " x " << n << endl;
        auto B = getLaplacian2D(n);
        double phoB = cos(M_PI  / (n + 1));
        double omega = 2.0 / (1.0 + sqrt(1 - pow(phoB, 4)));
        cout << "Omega: " << omega << endl;

        cout << "\e[1;32mPreparing\e[0m preconditioner..." << endl;
        auto Q = SSORPreconditioner(B, omega);

        cout << "\e[1;32mRunning\e[0m PCG Iteration..." << endl;
        double h = 1.0 / (n + 1);
        ColVector b2 = discreteF(n) * (h * h);
        ColVector x2 = PCG(B, b2, zeros(n * n, 1), Q, eps);
        cout << "Residule: " << norm(B * x2 - b2) << endl;
        pcolor("heatmap.png", x2.reshape(n, n));
        cout << "------------------------------------------------------------" << endl;
    } 
    else {
        if(testID == 3)
            cout << "--------- MG-Preconditioned CG For 2D Laplacian -----------" << endl;
        else
            cout << "--------- AMG-Preconditioned CG For 2D Laplacian ----------" << endl;
        cout << "Grid size: " << n << " x " << n << endl;
        auto B = getLaplacian2D(n);

        cout << "\e[1;32mPreparing\e[0m preconditioner..." << endl;
        Preconditioner* Q = (testID == 3) ? 
                            (Preconditioner*) new Lap2DMGPreconditioner(n) : 
                            (Preconditioner*) new AMGPreconditioner(B);
        
        int timest = clock();
        cout << "\e[1;32mRunning\e[0m PCG Iteration..." << endl;
        double h = 1.0 / (n + 1);
        ColVector b2 = discreteF(n) * (h * h);
        ColVector x2 = PCG(B, b2, zeros(n * n, 1), *Q, eps);
        cout << "Residule: " << norm(B * x2 - b2) << endl;
        cout << "PCG running time: " << std::setprecision(3) << (double)(clock()-timest)/CLOCKS_PER_SEC << "s" << endl;
        cout << std::setprecision(6);
        pcolor("heatmap.png", x2.reshape(n, n));
        cout << "-----------------------------------------------------------" << endl;
    }
    return 0;
}