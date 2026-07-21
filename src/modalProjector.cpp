#include "ModalProjector.h"

void ModalProjector::project(const Geometry& leftGeom, const ModeData& leftModes,
                             const Geometry& rightGeom, const ModeData& rightModes,
                             const SurfaceLoad& load,
                             std::vector<double>& forceL,
                             std::vector<double>& forceR) const {
    forceL.assign(leftModes.nModes, 0.0);
    forceR.assign(rightModes.nModes, 0.0);
    for (int modeIndex = 0; modeIndex < leftModes.nModes; ++modeIndex)
        for (int i = 0; i < leftGeom.nsurfl; ++i) for (int j = 0; j < leftGeom.nsurfz; ++j) {
            const int id = leftGeom.surfp[i][j]; if (id < 0) continue;
            const auto& u = leftModes.modes[modeIndex][id];
            forceL[modeIndex] += load.fxL[i][j]*u.ux + load.fyL[i][j]*u.uy + load.fzL[i][j]*u.uz;
        }
    for (int modeIndex = 0; modeIndex < rightModes.nModes; ++modeIndex)
        for (int i = 0; i < rightGeom.nsurfl; ++i) for (int j = 0; j < rightGeom.nsurfz; ++j) {
            const int id = rightGeom.surfp[i][j]; if (id < 0) continue;
            const auto& u = rightModes.modes[modeIndex][id];
            forceR[modeIndex] += load.fxR[i][j]*u.ux + load.fyR[i][j]*u.uy + load.fzR[i][j]*u.uz;
        }
}
