#pragma once
#include <cmath>
#include <vector>

// Identical, deterministic frozen shapes for actual user-main and adapted code.
inline void legacyFixture(int kind, Geometry& g, State& st, ModeData& modes) {
    g.nsurfl=g.nxsup=7; g.nsurfz=5; g.nPoints=35; g.zmax=17.; g.xsup=12.;
    g.surfp.assign(7,std::vector<int>(5)); g.ymid.assign(5,0.);
    g.points.resize(35); st.disp.resize(35);
    modes.nModes=2; modes.nPoints=35;
    modes.modes.assign(2,std::vector<Displacement>(35));
    const double base[7]={2.,1.5,1.,.4,.6,1.,1.5};
    for (int i=0;i<7;++i) for (int j=0;j<5;++j) {
        const int id=i*5+j; g.surfp[i][j]=id;
        double h=kind==0?.6:base[i];
        if (kind==2) h+=.8*std::cos(j*M_PI/2.);
        if (kind==4) h+=.12*std::sin(.6*i+.7*j);
        if (kind==5) h=0.;
        const double x=i*2., z=j*17./4.;
        g.points[id]={x,-.5*h,z};
        st.disp[id].ux=x+(kind==3?.05*i*std::sin(i):0.);
        st.disp[id].uy=-.5*h; st.disp[id].uz=z;
        modes.modes[0][id].uy=1.;
        modes.modes[1][id].ux=.3; modes.modes[1][id].uy=-.1; modes.modes[1][id].uz=.2;
    }
}
