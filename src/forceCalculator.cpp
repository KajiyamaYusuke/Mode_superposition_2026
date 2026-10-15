#include "ForceCalculator.h"
#include <iostream>

ForceCalculator::ForceCalculator(const Geometry& geomL_, const Geometry& geomR_,
                                 const ModeData& mdL_, const ModeData& mdR_,
                                 State& stL_, State& stR_, const SimulationParams& sp_)
    : geomL(geomL_), geomR(geomR_),
      modeDataL(mdL_), modeDataR(mdR_),
      stateL(stL_), stateR(stR_), sp(sp_),
      fxL(surfaceLoad.fxL), fyL(surfaceLoad.fyL), fzL(surfaceLoad.fzL),
      fxR(surfaceLoad.fxR), fyR(surfaceLoad.fyR), fzR(surfaceLoad.fzR),
      harea(sections.area) {}

void ForceCalculator::initialize() {
    // --- 左右それぞれのサイズを取得 ---
    int nModes_L = modeDataL.nModes;

    int nModes_R = modeDataR.nModes;

    // --- 左声帯 (L) の変数初期化 ---
    surfaceLoad.resize(geomL, geomR);
    fiL.assign(nModes_L, 0.0);

    // --- 右声帯 (R) の変数初期化 ---
    fiR.assign(nModes_R, 0.0);

    // --- 内部バッファ (流路断面数は左右で共通と仮定) ---
    nxsup = geomL.nxsup;

    fdisXL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));
    fdisYL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));

    fdisXR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));
    fdisYR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));

    prevFdisXL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));
    prevFdisYL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));

    prevFdisXR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));
    prevFdisYR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));

    // --- 共有変数の初期化 ---
    constexpr int flowSectionCount = 50;
    sections.resize(geomL, geomR, flowSectionCount);
    psurf.assign(flowSectionCount, 0.0);
    Ug.assign(stateL.nSteps, 0.0);       // ステップ数はL/R共通

    currentUg = 0.0;
    currentPg = 0.0;
    previousUg = 0.0;
    contactFlag = false;
    flowModel.initialize(sp, geomL, stateL.nSteps);

    rho = sp.rho ;
    mu  = sp.mu ;
    std::cout << "[ForceCalculator] initialized (Asymmetric): "
              << "nModesL=" << nModes_L << ", nModesR=" << nModes_R
              << ", nxsup=" << nxsup
              << ", L_sub=" << sp.L_sub
              << ", L_vt=" << sp.L_vt
              << ", separation_method=downstream_area_ratio"
              << std::endl;
    std::cout << "[ChannelSections] " << flowSectionCount
              << " fixed flow planes; retained initial constriction at i="
              << sections.referenceConstrictionI
              << " (sigma=" << sections.referenceConstrictionSigma << ")"
              << std::endl;
}

void ForceCalculator::projectLoadsToModes() {
    modalProjector.project(geomL, modeDataL, geomR, modeDataR,
                           surfaceLoad, fiL, fiR);
}

void ForceCalculator::updateChannelSections() {
    sectionBuilder.build(geomL, stateL, geomR, stateR, sections);
}
