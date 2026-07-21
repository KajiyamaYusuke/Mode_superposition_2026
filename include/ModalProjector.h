#pragma once

#include <vector>
#include "Geometry.h"
#include "ModeData.h"
#include "SurfaceLoad.h"

class ModalProjector {
public:
    void project(const Geometry& leftGeom, const ModeData& leftModes,
                 const Geometry& rightGeom, const ModeData& rightModes,
                 const SurfaceLoad& load,
                 std::vector<double>& forceL,
                 std::vector<double>& forceR) const;
};
