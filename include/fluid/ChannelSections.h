#pragma once

#include <vector>
#include "Geometry.h"
#include "State.h"

// Instantaneous geometry seen by the 1-D glottal-flow model.  It deliberately
// contains no flow history and no structural state.
struct ChannelSections {
    // Inactive legacy compatibility field; x planes are spaced in the
    // instantaneous whole-surface common x range, without a material anchor.
    std::vector<double> sigma;
    std::vector<double> x;       // common x = const section locations [mm]
    std::vector<double> area;    // open area of every section [mm^2]
    // Integral of the positive local gap cubed, ∫g_+^3 dz [mm^4].
    // It retains spanwise non-uniformity for the viscous flow resistance.
    std::vector<double> gapCubed;
    std::vector<bool> valid;     // geometric cut succeeded (not a physical closure)

    // Sign converting yRight-yLeft to the initially open-gap orientation.
    // The signed local gap is clamped to zero during area integration, so a
    // penetrated surface pair can never create a positive area by polygon
    // winding reversal.
    double gapSign = 1.0;

    // Inactive legacy anchor fields (reset by resize()).
    int referenceConstrictionI = -1;
    double referenceConstrictionSigma = 0.5;
    void resize(const Geometry& left, const Geometry& right, int nSections);
};

// Builds common flow sections by intersecting every constant-j surface line
// with a common x = const plane. Each side selects its inner intersection;
// only overlapping physical-z segments from consecutive valid j are used.
class ChannelSectionBuilder {
public:
    void build(const Geometry& leftGeom, const State& leftState,
               const Geometry& rightGeom, const State& rightState,
               ChannelSections& sections) const;
};
