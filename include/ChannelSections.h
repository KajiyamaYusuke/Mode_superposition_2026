#pragma once

#include <vector>
#include "Geometry.h"
#include "State.h"

// Instantaneous geometry seen by the 1-D glottal-flow model.  It deliberately
// contains no flow history and no structural state.
struct ChannelSections {
    // Fixed material-flow coordinate [0, 1].  These labels are selected once
    // from the undeformed mesh and are not re-sampled at later time steps.
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

    // The legacy i-grid cross-section with the smallest initial opening is
    // deliberately retained as one of sigma.  This prevents a coarse flow
    // grid from silently skipping an initially closed (or nearly closed)
    // constriction.
    int referenceConstrictionI = -1;
    double referenceConstrictionSigma = 0.5;
    void resize(const Geometry& left, const Geometry& right, int nSections);
};

// Builds common flow sections by intersecting every constant-j surface line
// with a common x = const plane.  The planes retain fixed material labels;
// their deformed x positions are updated once per structural state.
class ChannelSectionBuilder {
public:
    void build(const Geometry& leftGeom, const State& leftState,
               const Geometry& rightGeom, const State& rightState,
               ChannelSections& sections) const;
};
