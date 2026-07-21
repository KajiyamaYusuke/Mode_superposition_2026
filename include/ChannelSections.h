#pragma once

#include <vector>
#include "Geometry.h"
#include "State.h"

// Instantaneous geometry seen by the 1-D glottal-flow model.  It deliberately
// contains no flow history and no structural state.
struct ChannelSections {
    std::vector<double> sigma;   // fixed material-flow coordinate [0, 1]
    std::vector<double> x;       // common x = const section locations [mm]
    std::vector<double> area;    // open area of every section [mm^2]
    std::vector<bool> valid;     // geometric cut succeeded (not a physical closure)
    std::vector<std::vector<std::vector<double>>> degreeL;
    std::vector<std::vector<std::vector<double>>> degreeR;

    void resize(const Geometry& left, const Geometry& right, int nSections);
};

// Builds common flow sections by intersecting every constant-j surface line
// with a common x = const plane.  It must be called once per structural state.
class ChannelSectionBuilder {
public:
    void build(const Geometry& leftGeom, const State& leftState,
               const Geometry& rightGeom, const State& rightState,
               ChannelSections& sections) const;
};
