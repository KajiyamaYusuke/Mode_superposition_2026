#include "ModeData.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <regex>


void ModeData::initialize(int nModes_, const Geometry& geom) {
    nModes  = nModes_;
    int nPoints = geom.nPoints;

    // モード形状の初期化 [nModes][nPoints]
    modes.assign(nModes, std::vector<Displacement>(nPoints, Displacement()));

    // 固有振動数・減衰比の初期化
    frequencies.assign(nModes, 0.0);
    dampingRatios.assign(nModes, 0.0);
}

void ModeData::loadFromVTU(const std::string& filename, const Geometry& geom) {
    int nPoints = geom.nPoints;
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open VTU file: " + filename);

    std::string line;
    // 周波数ごとに {X, Y, Z} データをまとめる辞書
    std::map<std::string, std::vector<double>> xData, yData, zData;

    // X/Y/Z 成分を判定する正規表現
    std::regex rxX("Name=\"Displacement_field,_X-?component_@_([^\"]+)\"");
    std::regex rxY("Name=\"Displacement_field,_Y-?component_@_([^\"]+)\"");
    std::regex rxZ("Name=\"Displacement_field,_Z-?component_@_([^\"]+)\"");

    std::string currentKey;
    std::vector<double>* currentBuffer = nullptr;
    
    // 【追加】重複回避のためのマップ: 元の名前 -> 現在のサフィックス付き名前
    std::map<std::string, std::string> activeKeyMap; 

    while (std::getline(file, line)) {
        std::smatch match;

        // ... 前略 ...
if (std::regex_search(line, match, rxX)) {
    std::string rawKey = match[1]; // 例: "232.91_Hz"
    
    // --- デバッグ出力追加: 何が見つかったか表示 ---
    // std::cout << "[DEBUG] Found X-key: " << rawKey;

    // --- 重複チェックロジック ---
        if (xData.count(rawKey) && !xData[rawKey].empty()) {
            // 重複検出
            int counter = 2;
            std::string newKey = rawKey + "_v" + std::to_string(counter);
            while (xData.count(newKey) && !xData[newKey].empty()) {
                counter++;
                newKey = rawKey + "_v" + std::to_string(counter);
            }
            
            // --- デバッグ出力: 重複があったことを表示 ---
            std::cout << " -> COLLISION! Renaming to: " << newKey << std::endl;

            activeKeyMap[rawKey] = newKey; 
        } else {
            // --- デバッグ出力: 新規であることを表示 ---
            //std::cout << " -> New mode (Registered)" << std::endl;
            
            activeKeyMap[rawKey] = rawKey;
        }
        
        currentKey = activeKeyMap[rawKey];
        xData[currentKey] = {};
        currentBuffer = &xData[currentKey];
        continue;

            
        } else if (std::regex_search(line, match, rxY)) {
            std::string rawKey = match[1];
            // X成分で決定されたキー（サフィックス付きかもしれない）を使用
            if (activeKeyMap.find(rawKey) == activeKeyMap.end()) {
                 // XがなくYがいきなり来た場合の安全策（通常ありえない）
                 activeKeyMap[rawKey] = rawKey;
            }
            currentKey = activeKeyMap[rawKey];
            
            yData[currentKey] = {};
            currentBuffer = &yData[currentKey];
            continue;
            
        } else if (std::regex_search(line, match, rxZ)) {
            std::string rawKey = match[1];
            // X成分で決定されたキーを使用
            if (activeKeyMap.find(rawKey) == activeKeyMap.end()) {
                 activeKeyMap[rawKey] = rawKey;
            }
            currentKey = activeKeyMap[rawKey];

            zData[currentKey] = {};
            currentBuffer = &zData[currentKey];
            continue;
            
        } else if (line.find("</DataArray>") != std::string::npos) {
            currentBuffer = nullptr;
            continue;
        }

        // DataArray内部の数値読み込み
        if (currentBuffer) {
            std::istringstream ss(line);
            double v;
            while (ss >> v) currentBuffer->push_back(v);
        }
    }

    // --- モード名（周波数）の一覧を抽出 ---
    std::vector<std::string> modeKeys;
    for (auto& [key, _] : xData) modeKeys.push_back(key);

    std::sort(modeKeys.begin(), modeKeys.end(), [](const std::string& a, const std::string& b) {
        auto parseHz = [](const std::string& s) {
            // "_v2" などがついている場合に対応するため、"_Hz"を探してそこまでを数値化
            return std::stod(s.substr(0, s.find("_Hz"))); 
        };
        // 数値が同じ（重複モード）の場合は、名前に含まれる "_v" の有無などで順序を保つ
        double va = parseHz(a);
        double vb = parseHz(b);
        if (std::abs(va - vb) < 1e-9) return a < b; // 値が同じなら辞書順 ("_Hz" < "_Hz_v2")
        return va < vb;
    });

    nModes = modeKeys.size();
    modes.resize(nModes);

    for (int m = 0; m < nModes; ++m) {
        const std::string& key = modeKeys[m];

        const auto& x = xData[key];
        const auto& y = yData[key];
        const auto& z = zData[key];

        if (x.size() != nPoints || y.size() != nPoints || z.size() != nPoints) {
            std::cerr << "Error in Mode: " << key << " (Points expected: " << nPoints << ")" << std::endl;
            std::cerr << " Sizes -> X: " << x.size() << ", Y: " << y.size() << ", Z: " << z.size() << std::endl;
            throw std::runtime_error("Mode " + key + ": size mismatch with geometry points");
        }

        std::vector<Displacement> disp(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            disp[i] = { z[i], x[i], y[i] }; // 座標入れ替え
        }
        modes[m] = disp;
    }

    // std::cout << "\n=== LOADED MODE KEYS LIST (" << nModes << ") ===" << std::endl;
    // for (const auto& key : modeKeys) {
    //     std::cout << key << std::endl;
    // }
    std::cout << "=======================================\n" << std::endl;

    std::cout << "ModeData: Loaded " << nModes 
              << " modes (" << nPoints << " points each)." << std::endl;
}


void ModeData::loadFreqDamping(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open frequency file: " + filename);

    frequencies.clear();
    dampingRatios.clear();

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') continue;

        std::istringstream ss(line);
        double freqHz, omega, damping, Q;
        ss >> freqHz >> omega >> damping >> Q;

        frequencies.push_back(freqHz);
        dampingRatios.push_back(damping);
    }

    if (frequencies.size() != nModes) {
        throw std::runtime_error("Mode count mismatch between VTU and frequency file");
    }

    nModes = frequencies.size();

    std::cout << "ModeData: Loaded frequencies for " << nModes << " modes." << std::endl;
}

void ModeData::normalizeModes(double mass, const Geometry& geom){
    
    int nPoints = geom.nPoints;
    double nMass = mass / nPoints ;
    double ci;
    //std::cout<<"nMass = "<<nMass<<"\n";

    for (int imode=0; imode<nModes; imode++){
        ci = 0.0;
        for (int j=0; j<nPoints; j++){
            ci += modes[imode][j].ux*modes[imode][j].ux
                + modes[imode][j].uy*modes[imode][j].uy
                + modes[imode][j].uz*modes[imode][j].uz;
        }
        ci = 1.0 / std::sqrt( nMass * ci );
    
        for (int j=0; j<nPoints; j++){
            modes[imode][j].ux *= ci;
            modes[imode][j].uy *= ci;
            modes[imode][j].uz *= ci;
        }
    }
    
}
