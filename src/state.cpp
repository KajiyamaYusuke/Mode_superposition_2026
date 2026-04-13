#include "State.h"
#include <iostream>
#include <cmath>

void State::initialize(int nPoints_, int nModes_, int nSteps_, const Geometry& geom) {
    nPoints = nPoints_;
    nModes  = nModes_;
    nSteps  = nSteps_;

    // モード座標の初期化
    q.assign(nModes, 0.0);
    qdot.assign(nModes, 0.0);
    qddot.assign(nModes, 0.0);

    // 節点変位
    disp.assign(nPoints, Displacement());
    predictedDisp.assign(nPoints, Displacement());
    vel.assign(nPoints, Displacement());

    for(int i = 0; i < nPoints; i++){
        disp[i].ux = 0.0 + geom.points[i].x;
        disp[i].uy = 0.0 + geom.points[i].y;
        disp[i].uz = 0.0 + geom.points[i].z;
    }

    for(int i = 0; i < nPoints; i++){
        vel[i].ux = 0.0 ;
        vel[i].uy = 0.0 ;
        vel[i].uz = 0.0 ;
    }
}

void State::mode2uf(const Geometry& geom, const ModeData& modeData, int step) {
    if (step < 0 || step >= nSteps) return;

    for (int i = 0; i < nPoints; ++i) {
        predictedDisp[i].ux = 0.0;
        predictedDisp[i].uy = 0.0;
        predictedDisp[i].uz = 0.0;

        vel[i].ux = 0.0;
        vel[i].uy = 0.0;
        vel[i].uz = 0.0;

        //std::cout<<"after_mode2uf|disp[1]= "<<disp[1].ux<<std::endl;

        for (int m = 0; m < nModes; ++m) {
            double qi = qf[m];
            predictedDisp[i].ux += modeData.modes[m][i].ux * qi* 1.0e3;
            predictedDisp[i].uy += modeData.modes[m][i].uy * qi* 1.0e3;
            predictedDisp[i].uz += modeData.modes[m][i].uz * qi* 1.0e3;

            vel[i].ux += modeData.modes[m][i].ux * qfdot[m] * 1.0e3;
            vel[i].uy += modeData.modes[m][i].uy * qfdot[m]* 1.0e3;
            vel[i].uz += modeData.modes[m][i].uz * qfdot[m]* 1.0e3;
        }
        predictedDisp[i].ux += geom.points[i].x ;
        predictedDisp[i].uy += geom.points[i].y ;
        predictedDisp[i].uz += geom.points[i].z ;

    }

}





// state.cpp : uf2uの中
void State::uf2u() {

    q = qf;
    qdot = qfdot;
    qddot = qfddot;

    disp = predictedDisp; 
}

