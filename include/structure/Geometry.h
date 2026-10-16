#pragma once
#include <vector>
#include <array>
#include <string>
#include <map>

struct Point { double x=0.0, y=0.0, z=0.0; 
    Point operator-(const Point& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }
};

class Geometry {
public:
    int nPoints = 0;
    int nCells  = 0;
    int nsurfz  = 0;
    int nsurfl  = 0;
    int nxsup   = 0;    // surfArea で計算される補助値

    double zmax = 0.0;
    double xsup = 0.0;
    std::vector<double> ymid;

    std::vector<Point> points;
    std::vector<std::vector<int>> connect;
    std::vector<int> offsets;
    std::vector<int> types;
    std::array<int, 20> jtypes{}; // 配列サイズは適宜

    

    // surf関連
    std::map<int, Point> grid;           // GRID番号 → 座標
    std::vector<std::vector<int>> quads; // CQUAD4要素 → 節点番号リスト
    std::vector<std::vector<int>> surfp; // 出力: 格子点インデックス

    std::vector<int> surfl;
    std::vector<double> surflx, surfly, surflz;
    std::vector<std::vector<double>> sarea; 


    void loadFromVTK(const std::string& filename);
    void surfExtractFromNAS(const std::string& nasFile, int nsurfl_param, int nsurfz_param);
    void surfArea();

        // --- デバッグ / 出力 ---
    void print() const;
};
