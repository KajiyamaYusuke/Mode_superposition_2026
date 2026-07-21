
#pragma once
#include <array>
#include <fstream>
#include <utility>
#include <vector>
#include "Geometry.h"
#include "ModeData.h"
#include "State.h"
#include "SimulationParams.h"
#include "ChannelSections.h"
#include "SurfaceLoad.h"
#include "ModalProjector.h"
#include "FlowModel.h"

// Temporary façade for the coupled load assembly.  Geometry construction,
// acoustic state integration, load storage, and modal projection are owned by
// the focused components below; this class retains the established contact
// algorithm and public simulation-facing API.
class ForceCalculator {
private:
    // Declared before the compatibility references below so their storage is
    // constructed first.
    ChannelSections sections;
    ChannelSectionBuilder sectionBuilder;
    SurfaceLoad surfaceLoad;
    ModalProjector modalProjector;
    FlowModel flowModel;

public:
    ForceCalculator(const Geometry& geomL, const Geometry& geomR, 
                    const ModeData& mdL, const ModeData& mdR, 
                    State& stL, State& stR, const SimulationParams& sp);

    void initialize();
    // Coupling sequence called by Simulation.
    void updateChannelSections();
    void applyFluidLoads(double t, int n);
    void projectLoadsToModes();
    void applyContactLoads(int step = -1, int contactIter = -1);
    void setContactMonitor(int surfaceI, int surfaceJ);

    // Legacy placeholder; retained only because it is part of the old public
    // header.  Contact is calculated by applyContactLoads().
    void contactForce();
    void outputForceVectors(int step) const;
    void outputCorrespondenceOffsets(int step) const;

    std::vector<std::array<double,3>> traceL;
    std::vector<std::array<double,3>> traceR;

    // Compatibility-facing views.  Storage is owned by the focused components
    // below; retaining these names keeps the contact solver isolated during the
    // migration and avoids changing its numerics.
    std::vector<std::vector<double>>& fxL;
    std::vector<std::vector<double>>& fyL;
    std::vector<std::vector<double>>& fzL;
    std::vector<std::vector<double>> fdisL; // 左の接触力バッファ
    std::vector<double> fiL;                // 左のモード力
    std::vector<std::vector<std::vector<double>>>& degreeL;

    std::vector<std::vector<double>>& fxR;
    std::vector<std::vector<double>>& fyR;
    std::vector<std::vector<double>>& fzR;
    std::vector<std::vector<double>> fdisR; // 右の接触力バッファ
    std::vector<double> fiR;                // 右のモード力
    std::vector<std::vector<std::vector<double>>>& degreeR;

    std::vector<double> psurf;
    std::vector<double> Ug;        // Glottal flow history
    std::vector<double> minHarea;  // Minimum area 
    std::vector<double>& harea;  // 流路断面積

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
    // Deprecated compatibility implementation.  New calculations use
    // FlowModel through applyFluidLoads().
    void calcFlowStep(double t, double dt, double current_min_area);

public:

    // Ishizaka & Flanagan (1972) モデル用
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
    std::ofstream contactMonitorFile;
    std::ofstream contactIterationFile;
    std::ofstream contactSearchFile;
    std::ofstream flowDebugFile;
    int contactMonitorI = -1;
    int contactMonitorJ = -1;
    // Candidate pairs from the first contact iteration of the current time
    // step.  Reusing them in later iterations avoids repeated broad-phase
    // searches while retaining a full search at every new time step.
    std::vector<std::vector<std::pair<int, int>>> contactPairCache;

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
