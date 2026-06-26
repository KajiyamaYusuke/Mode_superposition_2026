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

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < nPoints; i++){
        disp[i].ux = 0.0 + geom.points[i].x;
        disp[i].uy = 0.0 + geom.points[i].y;
        disp[i].uz = 0.0 + geom.points[i].z;
    }

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < nPoints; i++){
        vel[i].ux = 0.0 ;
        vel[i].uy = 0.0 ;
        vel[i].uz = 0.0 ;
    }

    std::cout << "nxsup = " << geom.nxsup << std::endl;
    std::cout << "Last node Z coord = " << disp[geom.surfp[geom.nxsup - 1][0]].uz << std::endl;
}

void State::mode2uf(const Geometry& geom, const ModeData& modeData, int step) {
    if (step < 0 || step >= nSteps) return;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nPoints; ++i) {
        double predUx = 0.0;
        double predUy = 0.0;
        double predUz = 0.0;

        double velUx = 0.0;
        double velUy = 0.0;
        double velUz = 0.0;

        //std::cout<<"after_mode2uf|disp[1]= "<<disp[1].ux<<std::endl;

        for (int m = 0; m < nModes; ++m) {
            double qi = qf[m];
            double qdi = qfdot[m];
            const auto& modeAtPoint = modeData.modes[m][i];

            predUx += modeAtPoint.ux * qi* 1.0e3;
            predUy += modeAtPoint.uy * qi* 1.0e3;
            predUz += modeAtPoint.uz * qi* 1.0e3;

            velUx += modeAtPoint.ux * qdi * 1.0e3;
            velUy += modeAtPoint.uy * qdi* 1.0e3;
            velUz += modeAtPoint.uz * qdi* 1.0e3;
        }
        predictedDisp[i].ux = predUx + geom.points[i].x ;
        predictedDisp[i].uy = predUy + geom.points[i].y ;
        predictedDisp[i].uz = predUz + geom.points[i].z ;

        vel[i].ux = velUx;
        vel[i].uy = velUy;
        vel[i].uz = velUz;

    }

}





// state.cpp : uf2uの中
void State::uf2u() {

    q = qf;
    qdot = qfdot;
    qddot = qfddot;

    disp = predictedDisp; 
}
