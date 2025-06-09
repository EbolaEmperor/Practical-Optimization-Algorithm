/**
 * 改进版 GnuplotWrapper
 *   1. 构造时指定 PNG 文件名
 *   2. drawPoint / drawLine 只缓存数据
 *   3. 析构时一次性写出 plot 命令 + 多个 '-' 数据段
 *
 * 编译：g++ main.cpp -std=c++17 -o demo
 */

#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class GnuplotWrapper {
public:
    enum class Color { Red, Green, Blue, Magenta, Cyan, Black };

    explicit GnuplotWrapper(const std::string& filename,
                            int w = 1600, int h = 1200)
        : gp_(popen("gnuplot", "w"))
    {
        if (!gp_) throw std::runtime_error("无法启动 gnuplot");

        fprintf(gp_,
                "set terminal pngcairo size %d,%d enhanced\n"
                "set output '%s'\n"
                "set grid\n"
                "set size ratio -1\n",
                w, h, filename.c_str());
    }

    ~GnuplotWrapper() {
        if (!gp_) return;

        /* ---------- 1) 组织 plot 命令 ---------- */
        fprintf(gp_, "plot ");
        for (std::size_t i = 0; i < datasets_.size(); ++i) {
            const auto& d = datasets_[i];
            if (i) fputs(", ", gp_);
            fprintf(gp_, "'-' with %s %s %g lc rgb '%s' notitle",
                    d.is_line ? "lines"   : "points",
                    d.is_line ? "lw"      : "pt 7 ps",
                    d.thickness,
                    colorName_(d.col));
        }
        fputc('\n', gp_);

        /* ---------- 2) 发送所有数据 ---------- */
        for (const auto& d : datasets_) {
            for (auto [x, y] : d.pts)
                fprintf(gp_, "%f %f\n", x, y);
            fputs("e\n", gp_);
        }

        /* ---------- 3) 结束 ---------- */
        fputs("unset output\n", gp_);
        pclose(gp_);
    }

    void drawPoint(std::pair<double, double> p,
                   Color c = Color::Red, double size = 1.5)
    { datasets_.push_back({false, {p}, c, size}); }

    void drawLine(std::pair<double, double> p1,
                  std::pair<double, double> p2,
                  Color c = Color::Blue, double lw = 2.0)
    { datasets_.push_back({true, {p1, p2}, c, lw}); }

private:
    struct DataSet {
        bool is_line;                            // true=线, false=点
        std::vector<std::pair<double,double>> pts;
        Color col;
        double thickness;                        // 线宽或点大小
    };

    FILE* gp_{nullptr};
    std::vector<DataSet> datasets_;

    static const char* colorName_(Color c) {
        switch (c) {
            case Color::Red:     return "red";
            case Color::Green:   return "green";
            case Color::Blue:    return "blue";
            case Color::Magenta: return "magenta";
            case Color::Cyan:    return "cyan";
            default:             return "black";
        }
    }
};