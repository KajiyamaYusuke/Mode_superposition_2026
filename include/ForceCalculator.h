
#pragma once
#include <array>
#include <fstream>
#include <vector>
#include "Geometry.h"
#include "ModeData.h"
#include "State.h"
#include "SimulationParams.h"

class ForceCalculator {
public:
    ForceCalculator(const Geometry& geomL, const Geometry& geomR, 
                    const ModeData& mdL, const ModeData& mdR, 
                    State& stL, State& stR, const SimulationParams& sp);

    void initialize();
    void calcForce(double t, int n); // メイン計算関数
    
    // 追加: 接触力計算など既存の関数
    void contactForce();
    void calcDis(int step = -1, int contactIter = -1);
    void calcArea();
    void f2mode();
    void outputForceVectors(int step) const;
    void outputCorrespondenceOffsets(int step) const;

    std::vector<std::array<double,3>> traceL;
    std::vector<std::array<double,3>> traceR;

    std::vector<std::vector<double>> fxL, fyL, fzL;
    std::vector<std::vector<double>> fdisL; // 左の接触力バッファ
    std::vector<double> fiL;                // 左のモード力
    std::vector<std::vector<std::vector<double>>> degreeL;

    std::vector<std::vector<double>> fxR, fyR, fzR;
    std::vector<std::vector<double>> fdisR; // 右の接触力バッファ
    std::vector<double> fiR;                // 右のモード力
    std::vector<std::vector<std::vector<double>>> degreeR;

    std::vector<double> psurf;
    std::vector<double> Ug;        // Glottal flow history
    std::vector<double> minHarea;  // Minimum area 
    std::vector<double> harea;  // 流路断面積

    std::vector<double> Uu; // Upstream (Trachea) flow
    std::vector<double> Pu; // Upstream (Trachea) pressure
    std::vector<double> Ud; // Downstream (Vocal Tract) flow
    std::vector<double> Pd; // Downstream (Vocal Tract) pressure

    std::vector<std::vector<double>> prevFdisXL;
    std::vector<std::vector<double>> prevFdisYL;
    std::vector<std::vector<double>> prevFdisXR;
    std::vector<std::vector<double>> prevFdisYR;

    bool contactFlag;

    double currentUg;
    double currentPg;   // Subglottal pressure at glottis entry
    double max_force_diff;

    std::vector<std::vector<double>> contactForceL_ij;
    std::vector<std::vector<double>> contactForceR_ij;

    std::vector<std::array<double,3>> lineStartL;
    std::vector<std::array<double,3>> lineEndL;
    std::vector<std::vector<double>> fdisXL, fdisYL;
    std::vector<std::vector<double>> fdisXR, fdisYR;

    void resetPreviousContactForce();

private:
    // Ishizaka & Flanagan (1972) モデル用
    void calcFlowStep(double t, double dt, double current_min_area); 

    // 状態変数 (Previous step values)
    // Nsecg: subglottal sections, Nsecp: supraglottal sections

    
    // 現在のステップのUg
    double previousUg = 0.0;
    
    
    double currentPout; // Radiation pressure

    // 物理定数・回路定数 (Initializeで計算または設定)
    double rho, mu, c_sound; // 空気密度, 粘性係数, 音速
    double alpha1,alpha2, beta;
    double Lui, Cui, Lu, Cu, R2; // 気管系定数
    double La, Ca, Ra, Lr, Rr;   // 声道系定数
    double lg; // 声門長 (depth)

    bool hasVocalTract;
    std::ofstream debugForceFile;
    std::ofstream contactDebugFile;

    const Geometry& geomL;
    const Geometry& geomR;
    const ModeData& modeDataL;
    const ModeData& modeDataR;
    State& stateL;
    State& stateR;
    const SimulationParams& sp;
    int nxsup;

    double findMinHarea();
    int findNsep(double minH);
    int findNearestRightSurfacePointSameJ(int leftI, int j) const;

};
