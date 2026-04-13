#pragma once
#include <vector>
#include "ModeData.h"
#include "Geometry.h"
#include "Displacement.h"

class State {
public:
    int nPoints = 0;
    int nModes  = 0;
    int nSteps  = 0;   // タイムステップ数
    double dt   = 0.0; // 時間刻み幅

    std::vector<double> q, qdot, qddot; // [nModes][]
    std::vector<double> qf, qfdot, qfddot; 
    std::vector<Displacement> disp;                   // [nPoints]
    std::vector<Displacement> vel;
    std::vector<Displacement> predictedDisp;


    void initialize(int nPoints_, int nModes_, int nSteps, const Geometry& geom);
    void mode2uf(const Geometry& geom, const ModeData& modeData, int step);
    void uf2u();

    void calcArea(const Geometry& geom);
};

