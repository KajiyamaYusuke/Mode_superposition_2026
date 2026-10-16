#include "ModalProjector.h"

void ModalProjector::project(const Geometry& leftGeom, const ModeData& leftModes,
                             const Geometry& rightGeom, const ModeData& rightModes,
                             const SurfaceLoad& load,
                             std::vector<double>& forceL,
                             std::vector<double>& forceR) const {
    if (forceL.size() != static_cast<std::size_t>(leftModes.nModes))
        forceL.resize(leftModes.nModes);
    if (forceR.size() != static_cast<std::size_t>(rightModes.nModes))
        forceR.resize(rightModes.nModes);

    #pragma omp parallel
    {
        #pragma omp for nowait schedule(static)
        for (int modeIndex = 0; modeIndex < leftModes.nModes; ++modeIndex) {
            double sum = 0.0;
            for (int i = 0; i < leftGeom.nsurfl; ++i) for (int j = 0; j < leftGeom.nsurfz; ++j) {
                const int id = leftGeom.surfp[i][j]; if (id < 0) continue;
                const auto& u = leftModes.modes[modeIndex][id];
                sum += load.fxL[i][j]*u.ux + load.fyL[i][j]*u.uy + load.fzL[i][j]*u.uz;
            }
            forceL[modeIndex] = sum;
        }

        #pragma omp for schedule(static)
        for (int modeIndex = 0; modeIndex < rightModes.nModes; ++modeIndex) {
            double sum = 0.0;
            for (int i = 0; i < rightGeom.nsurfl; ++i) for (int j = 0; j < rightGeom.nsurfz; ++j) {
                const int id = rightGeom.surfp[i][j]; if (id < 0) continue;
                const auto& u = rightModes.modes[modeIndex][id];
                sum += load.fxR[i][j]*u.ux + load.fyR[i][j]*u.uy + load.fzR[i][j]*u.uz;
            }
            forceR[modeIndex] = sum;
        }
    }
}
