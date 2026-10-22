#include "LegacyGlottalMethods.h"
#include "FluidConstants.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

LegacyGlottal::PressureResult LegacyGlottal::wallPressure(const ChannelSections& s,
    const Geometry& geom, const GlottalMethodOptions& options, const SeparationSelection& selection,
    double rho, double mu, double ug, double pg, double rampedPs, double pd) {
    PressureResult out;
    out.pressure.assign(s.area.size(), 0.);
    auto& p=out.pressure;
    p[0]=pg;
    out.sepIndex=selection.sepIndex;
    // User main findNsep returns i+1, then load loop uses i<nsep-1.
    // Retain this deliberate exclusion of the constriction material row.
    int oldNsep=static_cast<int>(s.area.size());
    for (int i=1; i<static_cast<int>(s.area.size()); ++i)
        if (std::abs(s.area[i]-selection.minAreaMm2)<1e-8 || s.area[i]<=0.) { oldNsep=i+1; break; }
    out.loadEndExclusive=oldNsep-1;
    const bool open=selection.minAreaMm2>(options.separation?1e-6:kAreaClosedMm2);
    if (open) {
        const int end=options.separation?static_cast<int>(s.area.size())-1:selection.sepIndex;
        out.sepIndex=end;
        for (int i=1; i<=end; ++i) {
            if (!s.valid[i-1] || !s.valid[i]) { p[i]=p[i-1]; continue; }
            double dx=std::abs(s.x[i]-s.x[i-1]);
            if (options.area && options.wallViscosity) {
                // User-confirmed C++ uses the middle AP row of REFERENCE x.
                const int j=geom.nsurfz/2;
                dx=std::abs(geom.points.at(geom.surfp[i][j]).x-geom.points.at(geom.surfp[i-1][j]).x);
            }
            const double ap=std::max(s.area[i-1],1e-3);
            const double ac=std::max(s.area[i],1e-3);
            const double ugmm=ug*1e6;
            const double bernoulli=.5*rho*ugmm*ugmm*(1./(ap*ap)-1./(ac*ac));
            double viscous=0;
            if (options.wallViscosity) {
                const double h=(s.area[i]+s.area[i-1])/(2.*geom.zmax);
                viscous=12.*mu*dx/(geom.zmax*h*h*h)*ugmm*1e3;
            } else {
                const double g3=.5*(s.gapCubed[i]+s.gapCubed[i-1]);
                viscous=g3>1e-12?12.*mu*dx*ugmm/g3*1e3:0.;
            }
            p[i]=p[i-1]+bernoulli-viscous;
            if (!std::isfinite(p[i])) throw std::runtime_error("nonfinite switched wall pressure");
            if (p[i]>p[i-1]) {
                out.recoveryIndex=i;
                if (options.separation) {
                    // Source writes the first recovered value, then EXITs.
                    // Retain that value; all never-computed downstream rows
                    // remain zero when legacyPostSeparation is selected.
                    out.sepIndex=i;
                } else {
                    p[i]=p[i-1];
                    for (int k=i+1; k<=end; ++k) p[k]=p[i];
                }
                break;
            }
        }
    } else {
        if (options.separation) {
            out.sepIndex=oldNsep-1;
            for (int i=1; i<oldNsep-1; ++i) p[i]=pg;
            for (int i=oldNsep-1; i<static_cast<int>(p.size()); ++i) p[i]=pd;
        } else {
            for (int i=1; i<out.sepIndex; ++i) p[i]=pg;
            p[out.sepIndex]=pd;
        }
    }
    if (!options.postSeparation)
        for (int i=out.sepIndex+1; i<static_cast<int>(p.size()); ++i) p[i]=pd;
    return out;
}

double LegacyGlottal::interpolatePressure(const ChannelSections& s, const std::vector<double>& p,
    int sep, double x, double blend, double pd, bool zeroDownstream) {
    const double xs=s.x[sep];
    // A material-row plane is the mean of the mirrored x coordinates. At
    // an abrupt legacy zero-pressure boundary, a one-ulp L/R difference
    // must not assign p_sep on one fold and zero on the other. Resolve only
    // machine-equal coordinates to the upstream endpoint; this is not a
    // finite-length pressure blend and does not alter the modern pipeline.
    const double planeRoundoff=32.*std::numeric_limits<double>::epsilon()
        *std::max({1.,std::abs(x),std::abs(xs)});
    if (zeroDownstream && std::abs(x-xs)<=planeRoundoff) return p[sep];
    if (x<=xs) {
        if (x==xs) return p[sep];
        if (x<=s.x.front()) return p.front();
        const auto upper=std::upper_bound(s.x.begin(),s.x.begin()+sep+1,x);
        const int k=static_cast<int>(upper-s.x.begin())-1;
        if (k<0 || k>=sep || !(s.x[k+1]>s.x[k]))
            throw std::runtime_error("nonmonotonic material-row x cannot support modern physical-x pressure mapping");
        const double w=(x-s.x[k])/(s.x[k+1]-s.x[k]);
        return (1.-w)*p[k]+w*p[k+1];
    }
    if (zeroDownstream) return 0.;
    const double end=std::min(xs+blend,s.x.back());
    if (x<end && end>xs+1e-12) {
        const double t=std::clamp((x-xs)/(end-xs),0.,1.);
        const double smooth=t*t*(3.-2.*t);
        return (1.-smooth)*p[sep]+smooth*pd;
    }
    return pd;
}
