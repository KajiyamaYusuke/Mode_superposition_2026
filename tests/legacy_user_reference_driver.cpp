// Compiled only against archived, user-confirmed main headers and extracted
// VERBATIM methods by tools/kanaya_ablation/validate.py. No coupled solver run.
#include <bits/stdc++.h>
#define private public
#include "ForceCalculator.h"
#undef private
#include "legacy_snapshot_fixture.h"
#include "TimeIntegrator.h"

// Frozen airway state: original calcForce still executes its own pressure,
// findNsep and force code. Only time-advancing the airway is suppressed.
void ForceCalculator::calcFlowStep(double,double,double) {}

int main(int argc,char** argv) {
    if(argc!=2)return 2;
    std::ofstream out(argv[1]);
    out<<std::setprecision(17)<<"case,i,j,area,pressure,Fx,Fy,Fz,Q0,Q1,W\n";
    for(int kind=0;kind<6;++kind) {
        Geometry g; State st; ModeData md; legacyFixture(kind,g,st,md);
        st.harea.assign(7,0.); st.calcArea(g);
        SimulationParams params; params.iforce=0; params.rho=1.13; params.mu=1.83e-5;
        ForceCalculator f(g,md,st,params);
        f.nxsup=7; f.rho=params.rho; f.mu=params.mu; f.lg=17.;
        f.currentUg=3e-4; f.currentPg=1000.; f.Pd.assign(10,20.);
        f.Ug.assign(2,0.); f.minHarea.assign(2,0.); f.fi.assign(2,0.);
        f.psurf.assign(7,0.);
        f.fdis.assign(7,std::vector<double>(4,0.));
        f.calcForce(.1,0); f.f2mode();
        double tx=0,ty=0,tz=0;
        for(int i=0;i<7;++i)for(int j=0;j<5;++j) {
            double fx=f.fx[i][j],fy=f.fy[i][j],fz=f.fz[i][j],qi=fy,qj=.3*fx-.1*fy+.2*fz;
            out<<kind<<','<<i<<','<<j<<','<<st.harea[i]<<','<<f.psurf[i]<<','<<fx<<','<<fy<<','<<fz<<','
                <<qi<<','<<qj<<','<<qi*1e-6-qj*2e-6<<'\n';
            tx+=fx;ty+=fy;tz+=fz;
        }
        out<<kind<<",-1,-1,0,0,"<<tx<<','<<ty<<','<<tz<<','<<f.fi[0]<<','<<f.fi[1]<<','<<f.fi[0]*1e-6-f.fi[1]*2e-6<<'\n';
    }
    TimeIntegrator integrator;
    for(int k=0;k<16;++k) {
        double q,v,a;
        integrator.newmarkStep(k%2?0.:.2,std::sin(k)*1e-5,std::cos(k)*1e-3,k*.001,
            1e-5,(45.+k*10.)*2*M_PI,.06,.275625,.55,q,v,a);
        out<<100+k<<",-2,-2,0,0,"<<q<<','<<v<<','<<a<<",0,0,0\n";
    }
}
