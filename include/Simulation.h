
#include "SimulationParams.h"
#include "ForceCalculator.h"
#include "TimeIntegrator.h"
#include <iostream>

class Simulation {
public:
    SimulationParams params;
    Geometry geomL;
    ModeData mdataL;
    State    stateL;

    Geometry geomR;
    ModeData mdataR;
    State    stateR;

    ForceCalculator fCalc;
    TimeIntegrator integrator;
    fs::path runDir;

    Simulation();
    //     : fCalc(geom, mdata, state, params )            // 必須引数を渡して初期化// TimeIntegrator も同様
    // {}


    void initialize(const fs::path& parameterFile = "../input/param.txt");
    void run();
    void writeVTK(int step, const Geometry& geom, const State& state, const std::string& rdir, int nwrite);
    void writeVTKCombined(int step, const Geometry& geomL, const State& stateL, 
                      const Geometry& geomR, const State& stateR, 
                      const std::string& rdir, int nwrite);
    
};
