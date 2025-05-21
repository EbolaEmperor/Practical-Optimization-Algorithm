// 手搓的线性支持向量机（linear SVM）
// 速度不快，点太多算不过来

#include <bits/stdc++.h>
#include <random>
#include <fstream>
#define SILENCE
#include "lib/general_constraint.h"
using namespace std;

// C 大：更严格（可能过拟合）
// C 小：泛化性更好（可能欠拟合）
double C;
vector<tuple<double, double, int>> dtas;

// x = (w1, w2, b, xi1, xi2, ..., xiN)
double f(const ColVector &x){
    // 目标函数
    double res = 0.5 * (x[0]*x[0] + x[1]*x[1]);
    for (int i = 3; i < x.size(); i++)
        res += C * x[i];
    return res;
}

ColVector h(const ColVector &x){
    // 等式约束
    return ColVector();
}

ColVector g(const ColVector &x){
    // 不等式约束
    ColVector res(2 * dtas.size());
    ColVector w(2);
    w[0] = x[0];
    w[1] = x[1];
    for (int i = 0; i < dtas.size(); i++){
        double xi = x[0] * get<0>(dtas[i]) + x[1] * get<1>(dtas[i]) + x[2];
        res[i] = xi * get<2>(dtas[i]) + x[3 + i] - 1;
        res[dtas.size() + i] = x[3 + i];
    }
    return res;
}

int main(){
    int n;
    cout << "Please input the number of data points: " << endl;
    cin >> n;
    
    const double real_w1 = 1.14;
    const double real_w2 = 5.14;
    const double real_b = -1.919;
    const double error = 0.01;
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution_x(-0.2, 2.0);
    uniform_real_distribution<double> distribution_y(0.0, 0.4);
    normal_distribution<double> distribution_noise(0, error);

    ofstream file("data.txt");

    for(int i = 0; i < n; i++){
        double x = distribution_x(generator);
        double y = distribution_y(generator);
        double label = (x * real_w1 + y * real_w2 + real_b > 0) ? 1 : -1;
        // 添加噪声
        x += distribution_noise(generator);
        y += distribution_noise(generator) / 5;
        dtas.push_back({x, y, label});
        file << x << " " << y << " " << label << endl;
    }
    cout << "Please input the value of C: " << endl;
    cin >> C;

    // 初始化变量
    ColVector x0(2 + 1 + n);
    auto x = general_constraint_optimal_PHR(f, h, g, x0);

    // 输出结果
    cout << "Split hyperplane: " << endl;
    cout << x[0] << " * x + " << x[1] << " * y + " << x[2] << " = 0" << endl;
    file << x[0] << " " << x[1] << " " << x[2] << endl;
    return 0;
}