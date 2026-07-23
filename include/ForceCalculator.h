
#pragma once
#include <filesystem>
#include <fstream>
#include <limits>
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
public:
private:
    // Declared before the compatibility references below so their storage is
    // constructed first.
    ChannelSections sections;
    ChannelSectionBuilder sectionBuilder;
    SurfaceLoad surfaceLoad;
    ModalProjector modalProjector;
    FlowModel flowModel;
    std::filesystem::path outputDirectory_ = "../output";
    int lastSeparationIndex_ = -1;
    double lastSeparationX_ = std::numeric_limits<double>::quiet_NaN();
    double lastBlendEndX_ = std::numeric_limits<double>::quiet_NaN();
    double lastSeparationPressure_ = std::numeric_limits<double>::quiet_NaN();

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
    void setOutputDirectory(const std::filesystem::path& directory);
    const std::filesystem::path& outputDirectory() const { return outputDirectory_; }
    double sectionX(int index) const {
        return (index >= 0 && index < static_cast<int>(sections.x.size()))
            ? sections.x[index] : std::numeric_limits<double>::quiet_NaN();
    }
    double sectionGapCubed(int index) const {
        return (index >= 0 && index < static_cast<int>(sections.gapCubed.size()))
            ? sections.gapCubed[index] : std::numeric_limits<double>::quiet_NaN();
    }
    double outletPressure() const { return flowModel.outletPressure(); }
    int separationIndex() const { return lastSeparationIndex_; }
    double separationX() const { return lastSeparationX_; }
    double separationBlendEndX() const { return lastBlendEndX_; }
    double separationPressure() const { return lastSeparationPressure_; }

    // Compatibility-facing views.  Storage is owned by the focused components
    // below; retaining these names keeps the contact solver isolated during the
    // migration and avoids changing its numerics.
    std::vector<std::vector<double>>& fxL;
    std::vector<std::vector<double>>& fyL;
    std::vector<std::vector<double>>& fzL;
    std::vector<double> fiL;                // 左のモード力

    std::vector<std::vector<double>>& fxR;
    std::vector<std::vector<double>>& fyR;
    std::vector<std::vector<double>>& fzR;
    std::vector<double> fiR;                // 右のモード力

    std::vector<double> psurf;
    std::vector<double> Ug;        // Glottal flow history
    std::vector<double> minHarea;  // Minimum area 
    std::vector<double>& harea;  // 流路断面積

    std::vector<std::vector<double>> prevFdisXL;
    std::vector<std::vector<double>> prevFdisYL;
    std::vector<std::vector<double>> prevFdisXR;
    std::vector<std::vector<double>> prevFdisYR;

    bool contactFlag;

    double currentUg;
    double currentPg;   // Subglottal pressure at glottis entry
    double max_force_diff;
    double contact_force_residual = 0.0;
    double contact_penetration_residual = 0.0;
    double max_contact_penetration = 0.0;

    std::vector<std::vector<double>> contactForceL_ij;
    std::vector<std::vector<double>> contactForceR_ij;

    std::vector<std::vector<double>> fdisXL, fdisYL;
    std::vector<std::vector<double>> fdisXR, fdisYR;

    void resetPreviousContactForce();

    double previousUg = 0.0;
    // Local constants used only by wall-pressure reconstruction.  Acoustic
    // circuit state itself is owned exclusively by FlowModel.
    double rho = 0.0, mu = 0.0, lg = 0.0;
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
    double previous_contact_penetration_ = 0.0;

    double findMinHarea();
    int findNsep(double minH);
};
