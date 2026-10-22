#pragma once
#include "ChannelSections.h"
#include "GlottalMethodOptions.h"
#include "SeparationSelection.h"
#include "State.h"

namespace LegacyGlottal {
// Primary source: KajiyamaYusuke/Mode_superposition main (user-confirmed success).
// Public Yoshinaga Fortran is a secondary reference; differences are documented.
// Material-row comparisons are restricted to mirrored surfaces; no 50-plane resampling.
void buildArea(const Geometry&, const State&, const Geometry&, const State&, ChannelSections&);
SeparationSelection select(const ChannelSections&, const GlottalMethodOptions&, double areaRatio);
struct Coefficients { double Lg, Rv, Rk; };
Coefficients coefficients(double rho, double mu, double beta, double lisM, double lgM, double areaM2);
struct PressureResult {
    std::vector<double> pressure;
    int sepIndex = -1;
    int recoveryIndex = -1;
    int loadEndExclusive = -1;
};
PressureResult wallPressure(const ChannelSections&, const Geometry&, const GlottalMethodOptions&,
                            const SeparationSelection&, double rho, double mu, double flow,
                            double pg, double rampedPs, double pd);
Displacement areaVector(const Geometry&, const State&, bool isLeft, int i, int j, bool legacy);
double interpolatePressure(const ChannelSections&, const std::vector<double>&, int sepIndex,
                           double x, double blendMm, double pd, bool zeroDownstream);
}
