#include "LegacyGlottalMethods.h"
#include "FluidConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

void LegacyGlottal::buildArea(const Geometry& l, const State& ls, const Geometry& r,
                              const State& rs, ChannelSections& out) {
    if (l.nxsup != r.nxsup || l.nsurfz != r.nsurfz || l.nsurfl != r.nsurfl)
        throw std::runtime_error("legacy material-row area requires matching left/right topology");
    const int ni = l.nxsup, nj = l.nsurfz;
    if (ni < 3 || nj < 3) throw std::runtime_error("legacy surface too small");
    const double center = .5*(l.points.at(l.surfp[0][0]).y+r.points.at(r.surfp[0][0]).y);
    out.x.resize(ni); out.area.assign(ni, 0.); out.gapCubed.assign(ni, 0.);
    out.valid.assign(ni, true);
    for (int i = 0; i < ni; ++i) {
        out.x[i] = .5*(ls.disp.at(l.surfp[i][nj/2]).ux+rs.disp.at(r.surfp[i][nj/2]).ux);
        for (int j = 0; j < nj; ++j) {
            const auto& a = ls.disp.at(l.surfp[i][j]);
            const auto& b = rs.disp.at(r.surfp[i][j]);
            if (!std::isfinite(a.ux+a.uy+a.uz+b.ux+b.uy+b.uz)
                || std::abs(a.ux-b.ux)>1e-8 || std::abs(a.uz-b.uz)>1e-8
                || std::abs(a.uy+b.uy-2*center)>1e-8)
                throw std::runtime_error("legacy material-row method is restricted to mirror symmetry");
            if (j+1 == nj) continue;
            const auto& c = ls.disp.at(l.surfp[i][j+1]);
            const auto& d = rs.disp.at(r.surfp[i][j+1]);
            const double h0 = out.gapSign*(b.uy-a.uy), h1 = out.gapSign*(d.uy-c.uy);
            const double dz = std::abs(c.uz-a.uz);
            // User-confirmed C++ State::calcArea clamps each AP CELL's mean
            // half-gap, unlike public Fortran's signed trapezoid. Preserve it.
            out.area[i] += std::max(.5*(h0+h1),0.)*dz;
            double u=h0, v=h1, length=dz;
            if (u<=0 && v<=0) continue;
            if (u<0) { length*=v/(v-u); u=0; }
            else if (v<0) { length*=u/(u-v); v=0; }
            // A positive-gap integral supports an explicitly modern coefficient
            // comparison; it does not redefine the legacy cell-mean area.
            out.gapCubed[i] += .25*(u*u*u+u*u*v+u*v*v+v*v*v)*length;
        }
    }
}

SeparationSelection LegacyGlottal::select(const ChannelSections& s,
                                         const GlottalMethodOptions& options, double ratio) {
    if (!options.area)
        return selectSeparationByAreaRatio(s.area, s.valid, ratio, kAreaClosedMm2);
    SeparationSelection result;
    for (int i=0; i<static_cast<int>(s.area.size()); ++i) {
        if (!s.valid[i] || !std::isfinite(s.area[i])) continue;
        if (result.minIndex<0 || s.area[i]<result.minAreaMm2) {
            result.minIndex=i; result.minAreaMm2=s.area[i];
        }
    }
    if (result.minIndex<0) throw std::runtime_error("no legacy area");
    result.targetAreaMm2=ratio*result.minAreaMm2;
    result.sepIndex=result.minIndex;
    if (result.minAreaMm2<=kAreaClosedMm2) {
        for (int i=1; i<static_cast<int>(s.area.size()); ++i)
            if (s.valid[i] && s.area[i]<=0) { result.sepIndex=i; break; }
    } else {
        for (int i=result.minIndex; i<static_cast<int>(s.area.size()); ++i)
            if (s.valid[i] && s.area[i]>=result.targetAreaMm2) { result.sepIndex=i; break; }
    }
    result.sepAreaMm2=s.area[result.sepIndex];
    result.reason="legacy_material_row_minimum";
    return result;
}
