#include "LegacyGlottalMethods.h"
#include "ForceCalculator.h"
#ifdef ENABLE_KANAYA_DIAGNOSTICS
#include "KanayaDiagnostics.h"
#endif
#include <algorithm>
#include <cmath>
#include <stdexcept>

Displacement LegacyGlottal::areaVector(const Geometry& g, const State& st,
                                       bool left, int i, int j, bool legacy) {
    const auto& im=st.disp.at(g.surfp[i-1][j]); const auto& ip=st.disp.at(g.surfp[i+1][j]);
    const auto& jm=st.disp.at(g.surfp[i][j-1]); const auto& jp=st.disp.at(g.surfp[i][j+1]);
    const double dx=.5*(ip.ux-im.ux), dy=.5*(ip.uy-im.uy), iz=.5*(ip.uz-im.uz);
    const double jx=.5*(jp.ux-jm.ux), jy=.5*(jp.uy-jm.uy), dz=.5*(jp.uz-jm.uz);
    if (legacy) {
        // atan(dy/dx), not atan2; deterministic limits only for degenerate
        // synthetic edges. All production reference dx,dz are nonzero.
        const auto angle=[](double a,double b){return b!=0?std::atan(a/b):(a==0?0:std::copysign(M_PI/2,a));};
        const double tx=angle(dy,dx), tz=angle(jy,dz), ds=std::sqrt(dx*dx+dy*dy);
        const double sign=left?1.:-1.;
        return {sign*ds*dz*1e-6*std::cos(tz)*std::sin(tx),
                -sign*ds*dz*1e-6*std::cos(tz)*std::cos(tx),
                sign*ds*dz*1e-6*std::sin(tz)};
    }
    Displacement a{(dy*dz-iz*jy)*1e-6,(iz*jx-dx*dz)*1e-6,(dx*jy-dy*jx)*1e-6};
    if ((left && a.uy>0)||(!left && a.uy<0)) { a.ux=-a.ux; a.uy=-a.uy; a.uz=-a.uz; }
    return a;
}

void ForceCalculator::applySwitchedFluidLoads(double t,int n) {
    if (sp.iforce!=0) throw std::runtime_error("legacy fluid strategies require iforce=0");
    for (auto* field:{&fxL,&fyL,&fzL,&fxR,&fyR,&fzR})
        for (auto& row:*field) std::fill(row.begin(),row.end(),0.);
    lastSeparationSelection_=LegacyGlottal::select(sections,sp.glottalMethods,sp.flowSeparationAreaRatio);
    lastFlowTerminalIndex_=lastSeparationSelection_.sepIndex;
    previousUg=n>0?Ug[n-1]:0.;
    flowModel.advance(t,sp.dt,lastSeparationSelection_.minAreaMm2*1e-6,previousUg,sections);
    currentUg=flowModel.flowRate(); currentPg=flowModel.glottalPressure(); Ug[n]=currentUg;
    if (!std::isfinite(currentUg+currentPg)) throw std::runtime_error("nonfinite switched flow");
    const double pd=flowModel.downstreamPressure().front();
    const auto pressure=LegacyGlottal::wallPressure(sections,geomL,sp.glottalMethods,
        lastSeparationSelection_,rho,mu,currentUg,currentPg,flowModel.lungPressure(),pd);
    psurf=pressure.pressure;
    const int sep=pressure.sepIndex;
    lastPressureRecoveryIndex_=pressure.recoveryIndex;
    lastLegacyLoadEndExclusive_=pressure.loadEndExclusive;
    lastSeparationIndex_=sep; lastSeparationX_=sections.x.at(sep);
    lastBlendEndX_=sp.glottalMethods.postSeparation?lastSeparationX_:
        std::min(lastSeparationX_+sp.flowBlendLengthMm,sections.x.back());
    lastSeparationPressure_=psurf.at(sep);
    lastSeparationSelection_.sepIndex=sep;
    lastSeparationSelection_.sepAreaMm2=sections.area.at(sep);
    for (int side=0; side<2; ++side) {
        const auto& g=side==0?geomL:geomR; const auto& st=side==0?stateL:stateR;
        auto& fx=side==0?fxL:fxR; auto& fy=side==0?fyL:fyR; auto& fz=side==0?fzL:fzR;
        for (int i=1; i<g.nxsup-1; ++i) for (int j=1; j<g.nsurfz-1; ++j) {
            const int id=g.surfp[i][j];
            if (id<0 || g.surfp[i-1][j]<0 || g.surfp[i+1][j]<0 || g.surfp[i][j-1]<0 || g.surfp[i][j+1]<0)
                throw std::runtime_error("legacy comparison requires complete material-grid surface");
            double p=0.;
            if (sp.glottalMethods.surfacePressureLoad) {
                p=i<pressure.loadEndExclusive?psurf.at(i):0.;
                if (!sp.glottalMethods.postSeparation && st.disp[id].ux>lastSeparationX_)
                    p=i<pressure.loadEndExclusive?LegacyGlottal::interpolatePressure(sections,psurf,sep,st.disp[id].ux,sp.flowBlendLengthMm,pd,false):0.;
            } else p=LegacyGlottal::interpolatePressure(sections,psurf,sep,st.disp[id].ux,
                        sp.flowBlendLengthMm,pd,sp.glottalMethods.postSeparation);
            const auto a=LegacyGlottal::areaVector(g,st,side==0,i,j,sp.glottalMethods.surfacePressureLoad);
            fx[i][j]=p*a.ux; fy[i][j]=p*a.uy; fz[i][j]=p*a.uz;
#ifdef ENABLE_KANAYA_DIAGNOSTICS
            KanayaDiagnostics::node(side==0,i,j,p,a.ux,a.uy,a.uz,fx[i][j],fy[i][j],fz[i][j]);
#endif
        }
    }
}
