
#pragma once
#include <vector>
#include "Geometry.h"
#include "ModeData.h"
#include "State.h"
#include "SimulationParams.h"

class ForceCalculator {
public:
    ForceCalculator(const Geometry& geom, const ModeData& md, State& st, const SimulationParams& sp);

    void initialize();
    void calcForce(double t, int n); // メイン計算関数
    
    // 追加: 接触力計算など既存の関数
    void contactForce();
    void calcDis();
    void f2mode();
    void outputForceVectors(int step) const;

    // 変数へのアクセス用（必要に応じて）
    const std::vector<double>& getPsurf() const { return psurf; }

public:
    // --- 既存のメンバ ---
    std::vector<std::vector<double>> fx, fy, fz;
    std::vector<double> fi;
    std::vector<std::vector<double>> fdis; // dissipation force
    std::vector<double> psurf;
    std::vector<double> Ug;        // Glottal flow history
    std::vector<double> minHarea;  // Minimum area history
    bool contactFlag;

    double currentUg;

private:
    // --- 新規追加: 圧縮性流体モデル用の変数 ---
    // Ishizaka & Flanagan (1972) モデル用
    void calcFlowStep(double t, double dt, double current_min_area); // FortranのcalcFlow相当

    // 状態変数 (Previous step values)
    // Nsecg: subglottal sections, Nsecp: supraglottal sections
    static const int Nsecg = 3; // 気管のセクション数
    static const int Nsecp = 10; // 声道のセクション数（適宜調整）

    std::vector<double> Uu; // Upstream (Trachea) flow
    std::vector<double> Pu; // Upstream (Trachea) pressure
    std::vector<double> Ud; // Downstream (Vocal Tract) flow
    std::vector<double> Pd; // Downstream (Vocal Tract) pressure

    
    // 現在のステップのUg
    double previousUg = 0.0;
    
    double currentPg;   // Subglottal pressure at glottis entry
    double currentPout; // Radiation pressure

    // 物理定数・回路定数 (Initializeで計算または設定)
    double rho, mu, c_sound; // 空気密度, 粘性係数, 音速
    double alpha1,alpha2, beta;
    double Lui, Cui, Lu, Cu, R2; // 気管系定数
    double La, Ca, Ra, Lr, Rr;   // 声道系定数
    double lg; // 声門長 (depth)

    // --- 既存の参照 ---
    const Geometry& geom;
    const ModeData& modeData;
    State& state;
    const SimulationParams& sp;
    int nxsup;

    double findMinHarea();
    int findNsep(double minH);
};