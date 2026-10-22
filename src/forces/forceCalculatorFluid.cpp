#include "ForceCalculator.h"
#ifdef ENABLE_KANAYA_DIAGNOSTICS
#include "KanayaDiagnostics.h"
#endif
#include "FluidConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

void ForceCalculator::applyFluidLoads(double t, int n) {
    if (sp.glottalMethods.any()) {
        applySwitchedFluidLoads(t, n);
        return;
    }
    // 左右で共通の分割数を使用
    int nsurfz = geomL.nsurfz;
    const int nSurfaceSections = geomL.nxsup;
    const int nFlowSections = static_cast<int>(harea.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fxL.size()); ++i) {
        std::fill(fxL[i].begin(), fxL[i].end(), 0.0);
        std::fill(fyL[i].begin(), fyL[i].end(), 0.0);
        std::fill(fzL[i].begin(), fzL[i].end(), 0.0);
    }

    // 右側を右側サイズでクリア
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fxR.size()); ++i) {
        std::fill(fxR[i].begin(), fxR[i].end(), 0.0);
        std::fill(fyR[i].begin(), fyR[i].end(), 0.0);
        std::fill(fzR[i].begin(), fzR[i].end(), 0.0);
    }

    if (sp.iforce == 1) {
        // ==== sin波加振 ====

        const int ni = std::min(static_cast<int>(fxL.size()),
                                static_cast<int>(fxR.size()));
        const int nj = std::min(static_cast<int>(fxL[0].size()),
                                static_cast<int>(fxR[0].size()));

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < ni - 1; i++) {
            for (int j = 1; j < nj - 1; j++) {
                const double f = sp.famp * std::sin(2.0 * M_PI * sp.forcef * t);
                if (sp.forceDirection == 0) {
                    fxL[i][j] = f;
                    fxR[i][j] = f;
                } else {
                    // Equal-magnitude forces in opposite gap directions.
                    // This probes the medial/lateral (uy) compliance while
                    // preserving mirror symmetry of the loading.
                    fyL[i][j] = -f;
                    fyR[i][j] =  f;
                }
            }
        }
    } else if (sp.iforce == 0) {
        // ==== 1D flow model ====
        // Select before advancing: FlowModel uses this same shared function,
        // area ratio and closure threshold for both Lg/Rv integration limits.
        lastSeparationSelection_ = selectSeparationByAreaRatio(
            harea, sections.valid, sp.flowSeparationAreaRatio, kAreaClosedMm2);
        if (lastSeparationSelection_.sepIndex < 0) {
            throw std::runtime_error("fluid time_s=" + std::to_string(t)
                + " step=" + std::to_string(n)
                + " minimum_section=" + std::to_string(lastSeparationSelection_.minIndex)
                + ": " + lastSeparationSelection_.reason);
        }
        const int sepIndex = lastSeparationSelection_.sepIndex;
        const double minA = lastSeparationSelection_.minAreaMm2;

        previousUg = (n > 0 && n-1 < (int)Ug.size()) ? Ug[n-1] : 0.0;

        // Both geometric flow coefficients stop at this same separation.
        flowModel.advance(t, sp.dt, minA * 1e-6, previousUg, sections);
        currentUg = flowModel.flowRate();
        currentPg = flowModel.glottalPressure();

        if (!std::isfinite(currentUg) || !std::isfinite(currentPg)) {
            throw std::runtime_error("nonfinite flow/pressure at time_s="
                + std::to_string(t) + " section=" + std::to_string(sepIndex));
        }
        Ug[n] = currentUg;

        // Construct attached-flow pressure up to the physical separation plane.
        const auto& tractPressure = flowModel.downstreamPressure();
        const double downstreamPressure = tractPressure.empty() ? 0.0 : tractPressure.front();
        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = currentPg;

        if (minA > kAreaClosedMm2) {
            for (int i = 1; i <= sepIndex; ++i) {
                if (!sections.valid[i - 1] || !sections.valid[i]) {
                    // Legacy pressure hold across an undefined adjacent cut.
                    psurf[i] = psurf[i - 1];
                    continue;
                }
                double dx = std::abs(sections.x[i] - sections.x[i-1]);
                double h_prev = std::max(harea[i-1], 1e-3);
                double h_curr = std::max(harea[i], 1e-3);
                const double g3 = 0.5 * (sections.gapCubed[i] + sections.gapCubed[i - 1]);

                double Ugm = currentUg * 1e6;

                double bernoulli_term = 0.5 * rho * Ugm * Ugm * (1.0/(h_prev*h_prev) - 1.0/(h_curr*h_curr));
                double viscous_term = g3 > 1.0e-12
                    ? 12.0 * mu * dx * Ugm / g3 * 1e3 : 0.0;

                psurf[i] = psurf[i-1] + bernoulli_term - viscous_term;

                if (psurf[i] > psurf[i-1]) {
                    // Numerical pressure recovery before the selected
                    // constriction must not leave zero-filled stations behind.
                    // Keep the last attached-flow value through xSep; the
                    // downstream surface load uses the legacy pressure blend.
                    psurf[i] = psurf[i - 1];
                    for (int k = i + 1; k <= sepIndex; ++k) psurf[k] = psurf[i];
                    break;
                }

            }
        } else {
            for (int i = 1; i < sepIndex; ++i) psurf[i] = currentPg;
            psurf[sepIndex] = downstreamPressure;
        }

        const double separationPressure = psurf[sepIndex];
        if (!std::isfinite(separationPressure)) {
            throw std::runtime_error("nonfinite Psep at time_s="
                + std::to_string(t) + " section=" + std::to_string(sepIndex));
        }
        for (int i = sepIndex + 1; i < nFlowSections; ++i) psurf[i] = downstreamPressure;

        // Surface nodes receive pressure from the actual deformed x = const
        // flow planes, rather than from their nearest material-index column.
        // This keeps the area calculation and the pressure-to-surface mapping
        // on the same set of planes when the two folds are tilted differently.
        const double xSep = sections.x[sepIndex];
        const double xBlendEnd = std::min(xSep + sp.flowBlendLengthMm, sections.x.back());
        lastSeparationIndex_ = sepIndex;
        lastSeparationX_ = xSep;
        lastBlendEndX_ = xBlendEnd;
        lastSeparationPressure_ = separationPressure;
        auto pressureAtSurfaceX = [&](double xNode) {
            if (!std::isfinite(xNode)) {
                throw std::runtime_error("nonfinite surface x at time_s="
                    + std::to_string(t) + " section=" + std::to_string(sepIndex));
            }
            if (xNode <= xSep) {
                if (xNode == xSep) return separationPressure;
                if (xNode <= sections.x.front()) return psurf.front();
                const auto upper = std::upper_bound(
                    sections.x.begin(), sections.x.begin() + sepIndex + 1, xNode);
                const int k = static_cast<int>(upper - sections.x.begin()) - 1;
                const double w = (xNode - sections.x[k]) / (sections.x[k + 1] - sections.x[k]);
                return (1.0 - w) * psurf[k] + w * psurf[k + 1];
            }
            // Restore the legacy wall-load transition, not the constant
            // downstream psurf output. Never overwrite the acoustic Pd state.
            if (xNode < xBlendEnd && xBlendEnd > xSep + 1.0e-12) {
                const double s = std::clamp((xNode - xSep) / (xBlendEnd - xSep), 0.0, 1.0);
                const double smooth = s * s * (3.0 - 2.0 * s);
                return (1.0 - smooth) * separationPressure + smooth * downstreamPressure;
            }
            return downstreamPressure;
        };

        struct Vector3 { double x, y, z; };
        auto applyPressureForce = [&](const Geometry& geom, const State& state,
                                      bool isLeft, int i, int j, double pressure,
                                      double& fx, double& fy, double& fz) {
            const int ip = geom.surfp[i + 1][j], im = geom.surfp[i - 1][j];
            const int jp = geom.surfp[i][j + 1], jm = geom.surfp[i][j - 1];
            if (ip < 0 || im < 0 || jp < 0 || jm < 0 || !std::isfinite(pressure)) return;
            const Vector3 ti{
                0.5 * (state.disp[ip].ux - state.disp[im].ux),
                0.5 * (state.disp[ip].uy - state.disp[im].uy),
                0.5 * (state.disp[ip].uz - state.disp[im].uz)};
            const Vector3 tj{
                0.5 * (state.disp[jp].ux - state.disp[jm].ux),
                0.5 * (state.disp[jp].uy - state.disp[jm].uy),
                0.5 * (state.disp[jp].uz - state.disp[jm].uz)};
            Vector3 dA{
                (ti.y * tj.z - ti.z * tj.y) * 1.0e-6,
                (ti.z * tj.x - ti.x * tj.z) * 1.0e-6,
                (ti.x * tj.y - ti.y * tj.x) * 1.0e-6};
            if ((isLeft && dA.y > 0.0) || (!isLeft && dA.y < 0.0)) {
                dA.x = -dA.x; dA.y = -dA.y; dA.z = -dA.z;
            }
            fx = pressure * dA.x;
            fy = pressure * dA.y;
            fz = pressure * dA.z;
#ifdef ENABLE_KANAYA_DIAGNOSTICS
            KanayaDiagnostics::node(isLeft, i, j, pressure, dA.x, dA.y, dA.z, fx, fy, fz);
#endif
        };

        for (int i = 1; i < nSurfaceSections - 1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {

                int pid_ip1_L = geomL.surfp[i+1][j];
                int pid_im1_L = geomL.surfp[i-1][j];
                int pid_jp1_L = geomL.surfp[i][j+1];
                int pid_jm1_L = geomL.surfp[i][j-1];

                if (pid_ip1_L >= 0 && pid_im1_L >= 0 && pid_jp1_L >= 0 && pid_jm1_L >= 0) {
                    const int pid_L = geomL.surfp[i][j];
                    const double pressure = pid_L >= 0
                        ? pressureAtSurfaceX(stateL.disp[pid_L].ux) : 0.0;
                    applyPressureForce(geomL, stateL, true, i, j, pressure,
                                       fxL[i][j], fyL[i][j], fzL[i][j]);

                }

                int pid_ip1_R = geomR.surfp[i+1][j];
                int pid_im1_R = geomR.surfp[i-1][j];
                int pid_jp1_R = geomR.surfp[i][j+1];
                int pid_jm1_R = geomR.surfp[i][j-1];

                if (pid_ip1_R >= 0 && pid_im1_R >= 0 && pid_jp1_R >= 0 && pid_jm1_R >= 0) {
                    const int pid_R = geomR.surfp[i][j];
                    const double pressure = pid_R >= 0
                        ? pressureAtSurfaceX(stateR.disp[pid_R].ux) : 0.0;
                    applyPressureForce(geomR, stateR, false, i, j, pressure,
                                       fxR[i][j], fyR[i][j], fzR[i][j]);

                }

            }
        }

    }

}

double ForceCalculator::findMinHarea() const {
    const auto selection = selectSeparationByAreaRatio(
        harea, sections.valid, sp.flowSeparationAreaRatio, kAreaClosedMm2);
    if (selection.sepIndex < 0)
        throw std::runtime_error("channel minimum: " + selection.reason);
    return selection.minAreaMm2;
}
