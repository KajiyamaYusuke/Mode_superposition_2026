#pragma once
#include <vector>
#include <array>
#include "Geometry.h"
#include "Displacement.h"



class ModeData {
public:
    int nModes = 0;
    int nPoints = 0;

    std::vector<std::vector<Displacement>> modes; // [nModes][nPoints]
    std::vector<double> frequencies;
    std::vector<double> dampingRatios;

    void initialize(int nModes_, const Geometry& geom);
    // VTU から変位を読み込む
    void loadFromVTU(const std::string& filename, const Geometry& geom);

    // テキストから周波数・ダンピングを読み込む
    void loadFreqDamping(const std::string& filename);

    // 統合メソッド
    void loadAll(const std::string& vtuFile, const std::string& freqFile, const Geometry& geom) {
        loadFromVTU(vtuFile, geom);
        loadFreqDamping(freqFile);
    }
    void normalizeModes(double mass, const Geometry& geom);
};
