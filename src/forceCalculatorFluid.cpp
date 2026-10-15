#include "ForceCalculator.h"
#include "FluidConstants.h"
#include <algorithm>
#include <cmath>
#include <limits>

void ForceCalculator::applyFluidLoads(double t, int n) {
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
        const double minA = findMinHarea();

        previousUg = (n > 0 && n-1 < (int)Ug.size()) ? Ug[n-1] : 0.0;

        // 安全な面積を使って流体計算を進める
        flowModel.advance(t, sp.dt, minA * 1e-6, previousUg, sections);
        currentUg = flowModel.flowRate();
        currentPg = flowModel.glottalPressure();

        lastSeparationSelection_ = selectSeparationByAreaRatio(
            harea, sections.valid, sp.flowSeparationAreaRatio, kAreaClosedMm2);
        // Invalid geometry remains diagnosable and uses the historical last
        // interior pressure station as a safe fallback.
        const int selectedIndex = lastSeparationSelection_.sepIndex >= 0
            ? lastSeparationSelection_.sepIndex
            : std::max(0, nFlowSections - 2);

        Ug[n] = currentUg;

        // Construct the pressure only up to the physical separation plane.
        // Downstream loading is mapped with a smooth physical-x blend below.
        const int sepIndex = std::clamp(selectedIndex, 0, nFlowSections - 1);
        const auto& tractPressure = flowModel.downstreamPressure();
        const double downstreamPressure = tractPressure.empty() ? 0.0 : tractPressure.front();
        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = currentPg;

        if (minA > kAreaClosedMm2) {
            for (int i = 1; i <= sepIndex; ++i) {
                if (!sections.valid[i - 1] || !sections.valid[i]) {
                    // A geometric cut failure is not a physical closure.  Do
                    // not turn it into an artificial pressure singularity.
                    psurf[i] = psurf[i - 1];
                    continue;
                }
                // Pressure stations use their deformed common-plane positions.
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
                    // physical-x downstream blend handles the separation.
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
            if (!std::isfinite(xNode) || !std::isfinite(xSep)) {
                return 0.0;
            }

            // Upstream of separation, interpolate only between the common
            // physical section planes.  No surface material index enters the
            // decision, so the two tilted folds share one load cutoff.
            if (xNode <= xSep) {
            int bestSegment = -1;
            double bestDistance = std::numeric_limits<double>::infinity();
                for (int k = 0; k < sepIndex; ++k) {
                const double x0 = sections.x[k];
                const double x1 = sections.x[k + 1];
                if (!std::isfinite(x0) || !std::isfinite(x1)) continue;
                if ((xNode - x0) * (xNode - x1) > 0.0) continue;
                const double midpointDistance = std::abs(xNode - 0.5 * (x0 + x1));
                if (midpointDistance < bestDistance) {
                    bestDistance = midpointDistance;
                    bestSegment = k;
                }
            }

            if (bestSegment >= 0) {
                const double x0 = sections.x[bestSegment];
                const double x1 = sections.x[bestSegment + 1];
                const double dx = x1 - x0;
                const double w = std::abs(dx) > 1.0e-12
                    ? std::clamp((xNode - x0) / dx, 0.0, 1.0) : 0.5;
                return (1.0 - w) * psurf[bestSegment] + w * psurf[bestSegment + 1];
            }

                if (xNode < sections.x.front()) {
                    return psurf.front();
                }
                return separationPressure;
            }

            // A two-plane smoothstep transition suppresses a non-physical
            // force cliff after separation while tending to the tract value.
            if (std::isfinite(xBlendEnd) && xNode < xBlendEnd
                && xBlendEnd > xSep + 1.0e-12) {
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
    double minimum = std::numeric_limits<double>::infinity();
    for (int i = 1; i + 1 < static_cast<int>(harea.size()); ++i) {
        if (!sections.valid[i] || !std::isfinite(harea[i])) continue;
        minimum = std::min(minimum, harea[i]);
    }
    return std::isfinite(minimum) ? minimum : 0.0;
}
