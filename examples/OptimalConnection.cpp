/*
    最优连接问题：平面上给定 n 个点，如何以最小的代价将它们连接起来？

    显然连的线都是直的，但不一定全都是给定点之间直接相连。
    例如等边三角形的三个顶点，将三个顶点到三角形重心连起来，才是最优的。
*/

#include <lib/conjugate_gradient.h>
#include <lib/gnuplot_wrapper.h>
#include <lib/quassi_newton.h>
#include <vector>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <random>
#include <utility>
using namespace std;
using Point = pair<double, double>;

// Prim 算法求完全图的最小生成树（复杂度 O(n^2)，适用于稠密图）
class Prim {
public:
    Prim(int n) : n(n), graph(n, vector<double>(n, INF)), min_edge(n, INF), in_mst(n, false) {
        for (int i = 0; i < n; ++i) graph[i][i] = 0;
    }

    void add_edge(int u, int v, double weight) {
        graph[u][v] = graph[v][u] = weight;
    }

    double minimum_spanning_tree(bool record_edges = false) {
        min_edge[0] = 0;
        double total_weight = 0;
        vector<int> parent;
        if (record_edges) {
            mst_edges.clear();
            parent.resize(n, -1);
        }
        for (int i = 0; i < n; ++i) {
            int u = -1;
            for (int j = 0; j < n; ++j)
                if (!in_mst[j] && (u == -1 || min_edge[j] < min_edge[u]))
                    u = j;
            in_mst[u] = true;
            total_weight += min_edge[u];
            if (record_edges && parent[u] != -1)
                mst_edges.emplace_back(parent[u], u);
            for (int v = 0; v < n; ++v)
                if (graph[u][v] != INF && !in_mst[v])
                    if (min_edge[v] > graph[u][v]) {
                        min_edge[v] = graph[u][v];
                        if (record_edges) parent[v] = u;
                    }
        }
        return total_weight;
    }
    
    // 添加一个方法来获取记录的MST边
    const vector<pair<int, int>>& get_mst_edges() const {
        return mst_edges;
    }

private:
    static constexpr double INF = numeric_limits<double>::max();
    int n;
    vector<vector<double>> graph;
    vector<double> min_edge;
    vector<bool> in_mst;
    vector<pair<int, int>> mst_edges;
};

inline double dist(const Point& a, const Point& b) {
    return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
}

vector<Point> fixed_points;

// 将给定点和自由点合并为一个点集，建完全图，求最小生成树
double costFunc(const ColVector &free_points, bool output) {
    int n = fixed_points.size();
    int m = free_points.size() / 2;
    auto points = fixed_points;
    for (int i = 0; i < m; ++i) {
        points.emplace_back(free_points[2 * i], free_points[2 * i + 1]);
    }
    int total_points = points.size();
    Prim prim(total_points);
    for (int i = 0; i < total_points; ++i) {
        for (int j = i + 1; j < total_points; ++j) {
            double weight = dist(points[i], points[j]);
            prim.add_edge(i, j, weight);
        }
    }
    double ans = prim.minimum_spanning_tree(output);
    
    if (output) {
        GnuplotWrapper plotter("optimal_connection.png");
        auto edges = prim.get_mst_edges();
        for (const auto& edge : edges)
            plotter.drawLine(points[edge.first], points[edge.second]);
        for (int i = 0; i < n; ++i)
            plotter.drawPoint(points[i], GnuplotWrapper::Color::Black, 1.5);
    }

    return ans;
}

double costFunc(const ColVector &free_points) {
    return costFunc(free_points, false);
}

namespace FermatPoint {

    inline Point sub(const Point& a, const Point& b) {
        return {a.first - b.first, a.second - b.second};
    }
    inline double dot(const Point& a, const Point& b) {
        return a.first * b.first + a.second * b.second;
    }
    inline double norm(const Point& v){
        return std::sqrt(dot(v, v));
    }
    inline bool obtuse120(const Point& u, const Point& v){
        double du = norm(u), dv = norm(v);
        if (du == 0 || dv == 0) return true;      // 退化情况
        return dot(u, v) <= -0.5 * du * dv;
    }

    /**
     * @brief 计算三角形的第一费马点
     * @param A,B,C 三角形顶点（pair<double,double> 存储）
     * @return 费马点坐标；若存在角 ≥120° 则返回对应顶点
     */
    Point FermatPoint(const Point& A, const Point& B, const Point& C, double thr = 1e18)
    {
        // --------- 判断是否存在 ≥120° 的顶角 ---------
        Point AB = sub(B, A), AC = sub(C, A);
        Point BA = sub(A, B), BC = sub(C, B);
        Point CA = sub(A, C), CB = sub(B, C);

        if (obtuse120(AB, AC)) throw std::invalid_argument("Angle A ≥ 120°");
        if (obtuse120(BA, BC)) throw std::invalid_argument("Angle B ≥ 120°");
        if (obtuse120(CA, CB)) throw std::invalid_argument("Angle C ≥ 120°");
        if (dist(A, B) >= thr && dist(B, C) >= thr && dist(C, A) >= thr)
            throw std::invalid_argument("Distance between points is too large");

        // --------- 否则用正弦重心公式求费马点 ---------
        const double PI = std::acos(-1.0);
        const double sin120 = std::sqrt(3.0) / 2.0;   // sin(120°)

        // 三边长度
        double a = norm(BC);   // 对应角 A
        double b = norm(CA);   // 对应角 B
        double c = norm(AB);   // 对应角 C

        // 利用余弦定理求三角形内角（弧度）
        auto safe_acos = [](double x)
        {
            x = std::max(-1.0, std::min(1.0, x));
            return std::acos(x);
        };
        double alpha = safe_acos((b*b + c*c - a*a) / (2.0 * b * c));
        double beta  = safe_acos((c*c + a*a - b*b) / (2.0 * c * a));
        double gamma = safe_acos((a*a + b*b - c*c) / (2.0 * a * b));

        // 权重 w = 1 / sin(120° - 角)
        double wA = 1.0 / std::sin(2.0 * PI / 3.0 - alpha);
        double wB = 1.0 / std::sin(2.0 * PI / 3.0 - beta );
        double wC = 1.0 / std::sin(2.0 * PI / 3.0 - gamma);

        double sum = wA + wB + wC;

        Point P;
        P.first  = (wA * A.first  + wB * B.first  + wC * C.first ) / sum;
        P.second = (wA * A.second + wB * B.second + wC * C.second) / sum;
        return P;
    }

    vector<Point> getFermatPoints(const vector<Point>& points, double thr = 1e18) {
        vector<Point> fermat_points;
        int n = points.size();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                for (int k = j + 1; k < n; k++) {
                    try {
                        // 尝试计算费马点
                        fermat_points.push_back(FermatPoint(points[i], points[j], points[k], thr));
                    } catch (const std::invalid_argument&) {
                        // 如果存在 ≥120° 的角，则不添加该点
                    }
                }
        return fermat_points;
    }
}

int main() {
    const int n = 25; // 固定点数量

    mt19937 gen0(114);
    uniform_real_distribution<double> dis(0.0, 100.0);
    for (int i = 0; i < n; ++i) {
        fixed_points.emplace_back(dis(gen0), dis(gen0));
    }

    // 寻找所有费马点，距离阈值为 30.0，超过这个距离的点构成的三角形不考虑
    auto fermat_points = FermatPoint::getFermatPoints(fixed_points, 30.0);
    cout << "Fermat points found: " << fermat_points.size() << endl;

    cout << "--------------------- FR Nonliner CG Method ---------------------" << endl;
    ColVector min_opt;
    double min_cost = numeric_limits<double>::max();

    for (int T = 0; T < 200; T++) {
        const int m = n; // 自由点数量
        ColVector free_points(2 * m);
        shuffle(fermat_points.begin(), fermat_points.end(), gen0);
        for (int i = 0; i < m; ++i) {
            free_points[2 * i]     = fermat_points[i].first;
            free_points[2 * i + 1] = fermat_points[i].second;
        }
        auto opt = CG_FR_gradfree(costFunc, free_points, 1e-3, 0.1, 0.4);
        double cost = costFunc(opt);
        if (cost < min_cost) {
            min_cost = cost;
            min_opt = opt;
        }
        cout << "Minimum cost: " << cost << " (Min: " << min_cost << ")" << endl;
    }
    min_cost = costFunc(min_opt, true);
    return 0;
}