#pragma once

#include <vector>
#include "Geometry.h"

// Nodal loads on one or both medial surfaces.  This is the only mutable load
// container shared by pressure and contact calculations.
struct SurfaceLoad {
    std::vector<std::vector<double>> fxL, fyL, fzL;
    std::vector<std::vector<double>> fxR, fyR, fzR;

    void resize(const Geometry& left, const Geometry& right);
    void clear();
};
