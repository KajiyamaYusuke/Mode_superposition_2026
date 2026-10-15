#include "ForceCalculator.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>

#ifndef PROFILE_CALCDIS_EVERY
#define PROFILE_CALCDIS_EVERY 1000
#endif

void ForceCalculator::applyContactLoads(int step, int contactIter, bool diagnosticOnly) {

    const bool savedContactFlag = contactFlag;
    const double savedMaxForceDiff = max_force_diff;
    const double savedForceResidual = contact_force_residual;
    const double savedPenResidual = contact_penetration_residual;
    const double savedMaxPen = max_contact_penetration;
    const double savedPreviousPen = previous_contact_penetration_;

#ifdef PROFILE_CALCDIS
    using Clock = std::chrono::high_resolution_clock;

    auto now = []() {
        return Clock::now();
    };

    auto elapsed_ms = [](const auto& t0, const auto& t1) {
        return std::chrono::duration<double, std::milli>(t1 - t0).count();
    };

    static long long prof_calls = 0;

    static double prof_clear_ms  = 0.0;
    static double prof_build_ms  = 0.0;
    static double prof_active_ms = 0.0;
    static double prof_pair_ms   = 0.0;
    static double prof_accum_ms  = 0.0;
    static double prof_total_ms  = 0.0;

    static long long prof_total_segL = 0;
    static long long prof_total_segR = 0;
    static long long prof_total_process_pairs = 0;
    static long long prof_total_contacts = 0;
    static long long prof_total_active_candidates = 0;
    static long long prof_max_activeR = 0;

    long long call_segL = 0;
    long long call_segR = 0;
    long long call_process_pairs = 0;
    long long call_contacts = 0;
    long long call_active_candidates = 0;
    long long call_max_activeR = 0;

    const auto prof_t_total0 = now();
    auto prof_t0 = now();
#endif

    contactFlag = false;
    max_force_diff = 0.0;

    double maxPenetration = 0.0;

    double apZMin = std::numeric_limits<double>::infinity();
    double apZMax = -std::numeric_limits<double>::infinity();
    for (const auto& point : geomL.points) {
        apZMin = std::min(apZMin, point.z);
        apZMax = std::max(apZMax, point.z);
    }
    auto contactIsExcluded = [&](int j, double zMm) {
        if (!sp.apContactDiagnosticEnabled) return false;
        if (sp.apContactExcludedDistanceMm > 0.0) {
            return zMm <= apZMin + sp.apContactExcludedDistanceMm
                || zMm >= apZMax - sp.apContactExcludedDistanceMm;
        }
        return sp.apContactExcludedRowsPerEnd > 0
            && (j < sp.apContactExcludedRowsPerEnd
                || j >= std::max(0, geomL.nsurfz - sp.apContactExcludedRowsPerEnd));
    };

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fdisXL.size()); ++i) {
        std::fill(fdisXL[i].begin(), fdisXL[i].end(), 0.0);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fdisYL.size()); ++i) {
        std::fill(fdisYL[i].begin(), fdisYL[i].end(), 0.0);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fdisXR.size()); ++i) {
        std::fill(fdisXR[i].begin(), fdisXR[i].end(), 0.0);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fdisYR.size()); ++i) {
        std::fill(fdisYR[i].begin(), fdisYR[i].end(), 0.0);
    }

#ifdef PROFILE_CALCDIS
    {
        const auto prof_t1 = now();
        prof_clear_ms += elapsed_ms(prof_t0, prof_t1);
    }
#endif

    // ------------------------------------------------------------
    // 1. Contact stiffness reference
    // ------------------------------------------------------------
    // The existing kc1/kc2 values are reduced legacy coefficients.  Their
    // conversion is deliberately tied to one explicit, fixed reference
    // frequency, never to the first participating mode.  Frequency sweeps
    // therefore change the structure only, not the contact law.
    const double contactOmega = 2.0 * M_PI * sp.contactReferenceFrequencyHz;
    const double contactOmega2 = contactOmega * contactOmega;

    constexpr double EPS_X   = 1.0e-12;
    constexpr double EPS_LEN = 1.0e-12;

    struct Pos {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    auto clamp01 = [](double v) {
        return std::max(0.0, std::min(1.0, v));
    };

    auto lerp = [](double a, double b, double s) {
        return (1.0 - s) * a + s * b;
    };

    auto interpPos = [&](const auto& disp, int pid0, int pid1, double s) -> Pos {
        return Pos{
            lerp(disp[pid0].ux, disp[pid1].ux, s),
            lerp(disp[pid0].uy, disp[pid1].uy, s),
            lerp(disp[pid0].uz, disp[pid1].uz, s)
        };
    };

    auto interpVelocity = [&](const auto& velocity, int pid0, int pid1, double s) -> Pos {
        return Pos{
            lerp(velocity[pid0].ux, velocity[pid1].ux, s),
            lerp(velocity[pid0].uy, velocity[pid1].uy, s),
            lerp(velocity[pid0].uz, velocity[pid1].uz, s)
        };
    };

    auto paramAtX = [&](double x, double x0, double x1) {
        const double den = x1 - x0;
        if (std::abs(den) < EPS_X) return 0.5;
        return clamp01((x - x0) / den);
    };

    auto norm2 = [](double x, double y) {
        return std::sqrt(x * x + y * y);
    };

    auto evalPredAtX =
        [&](int pidL0, int pidL1,
            int pidR0, int pidR1,
            double x,
            double& sL,
            double& sR,
            Pos& pL,
            Pos& pR) {
            const double xL0 = stateL.predictedDisp[pidL0].ux;
            const double xL1 = stateL.predictedDisp[pidL1].ux;
            const double xR0 = stateR.predictedDisp[pidR0].ux;
            const double xR1 = stateR.predictedDisp[pidR1].ux;

            sL = paramAtX(x, xL0, xL1);
            sR = paramAtX(x, xR0, xR1);

            pL = interpPos(stateL.predictedDisp, pidL0, pidL1, sL);
            pR = interpPos(stateR.predictedDisp, pidR0, pidR1, sR);

            return pR.y - pL.y;
        };

    // ------------------------------------------------------------
    // 2. Negative gap interval and centroid
    // ------------------------------------------------------------
    auto computePenetrationCentroid =
        [&](double xa, double xb,
            double gapA,
            double gapB,
            double& xc0,
            double& xc1,
            double& xContact,
            double& penMean,
            double& penMax) -> bool {
            if (xb <= xa) return false;

            if (gapA >= 0.0 && gapB >= 0.0) {
                return false;
            }

            xc0 = xa;
            xc1 = xb;

            if (gapA >= 0.0 && gapB < 0.0) {
                const double r = gapA / (gapA - gapB);
                xc0 = xa + r * (xb - xa);
                xc1 = xb;
            } else if (gapA < 0.0 && gapB >= 0.0) {
                const double r = gapA / (gapA - gapB);
                xc0 = xa;
                xc1 = xa + r * (xb - xa);
            } else {
                xc0 = xa;
                xc1 = xb;
            }

            if (xc1 <= xc0) return false;

            auto gapLinear = [&](double x) {
                const double u = (x - xa) / (xb - xa);
                return lerp(gapA, gapB, u);
            };

            const double p0 = std::max(0.0, -gapLinear(xc0));
            const double p1 = std::max(0.0, -gapLinear(xc1));

            const double L = xc1 - xc0;
            const double areaP = 0.5 * (p0 + p1) * L;

            if (areaP <= 0.0) return false;

            const double momentLocal =
                L * L * (p0 / 2.0 + (p1 - p0) / 3.0);

            xContact = xc0 + momentLocal / areaP;

            penMean = areaP / L;
            penMax  = std::max(p0, p1);

            return true;
        };

    // ------------------------------------------------------------
    // 3. XY contact normal
    // ------------------------------------------------------------
    auto computeXYNormal =
        [&](int pidL0, int pidL1,
            int pidR0, int pidR1,
            const Pos& pL,
            const Pos& pR,
            double& nx,
            double& ny) {
            const auto& L0 = stateL.predictedDisp[pidL0];
            const auto& L1 = stateL.predictedDisp[pidL1];
            const auto& R0 = stateR.predictedDisp[pidR0];
            const auto& R1 = stateR.predictedDisp[pidR1];

            double txL = L1.ux - L0.ux;
            double tyL = L1.uy - L0.uy;
            double lenL = norm2(txL, tyL);

            if (lenL > EPS_LEN) {
                txL /= lenL;
                tyL /= lenL;
            } else {
                txL = 0.0;
                tyL = 0.0;
            }

            double txR = R1.ux - R0.ux;
            double tyR = R1.uy - R0.uy;
            double lenR = norm2(txR, tyR);

            if (lenR > EPS_LEN) {
                txR /= lenR;
                tyR /= lenR;
            } else {
                txR = 0.0;
                tyR = 0.0;
            }

            double tx = txL + txR;
            double ty = tyL + tyR;
            double lenT = norm2(tx, ty);

            if (lenT <= EPS_LEN) {
                tx = txL;
                ty = tyL;
                lenT = norm2(tx, ty);
            }

            if (lenT <= EPS_LEN) {
                nx = 0.0;
                ny = 1.0;
                return;
            }

            tx /= lenT;
            ty /= lenT;

            nx = -ty;
            ny =  tx;

            // The normal orientation must be fixed by the *reference*
            // left-to-right separation.  Orienting it by pR-pL at the current
            // point would make (pR-pL)·n positive by construction and hence
            // make normal penetration mathematically impossible to detect.
            const double xL0m = 0.5 * (geomL.points[pidL0].x + geomL.points[pidL1].x);
            const double yL0m = 0.5 * (geomL.points[pidL0].y + geomL.points[pidL1].y);
            const double xR0m = 0.5 * (geomR.points[pidR0].x + geomR.points[pidR1].x);
            const double yR0m = 0.5 * (geomR.points[pidR0].y + geomR.points[pidR1].y);
            const double dxLR0 = xR0m - xL0m;
            const double dyLR0 = yR0m - yL0m;

            double lenLR0 = norm2(dxLR0, dyLR0);

            if (lenLR0 > EPS_LEN) {
                const double dot = nx * dxLR0 + ny * dyLR0;
                if (dot < 0.0) {
                    nx = -nx;
                    ny = -ny;
                }
            } else {
                if (ny < 0.0) {
                    nx = -nx;
                    ny = -ny;
                }
            }

            const double nlen = norm2(nx, ny);
            if (nlen > EPS_LEN) {
                nx /= nlen;
                ny /= nlen;
            } else {
                nx = 0.0;
                ny = 1.0;
            }
        };

    auto addSegmentForceXY =
        [&](std::vector<std::vector<double>>& fbufX,
            std::vector<std::vector<double>>& fbufY,
            int i0,
            int i1,
            int j,
            double s,
            double Fx,
            double Fy) {
            const double w0 = 1.0 - s;
            const double w1 = s;

            fbufX[i0][j] += w0 * Fx;
            fbufX[i1][j] += w1 * Fx;

            fbufY[i0][j] += w0 * Fy;
            fbufY[i1][j] += w1 * Fy;

        };

    // ------------------------------------------------------------
    // 4. Loop over same j and x-overlapping segment pairs
    // ------------------------------------------------------------
    const int nSegL = std::min({
        geomL.nxsup,
        static_cast<int>(fdisXL.size()),
        static_cast<int>(geomL.sarea.size())
    }) - 1;

    const int nSegR = std::min({
        geomR.nxsup,
        static_cast<int>(fdisXR.size()),
        static_cast<int>(geomR.sarea.size())
    }) - 1;

    const int maxJ = std::min({
        geomL.nsurfz,
        geomR.nsurfz,
        static_cast<int>(fdisXL.empty() ? 0 : fdisXL[0].size()),
        static_cast<int>(fdisXR.empty() ? 0 : fdisXR[0].size())
    });

    struct ContactSeg {
        int i    = -1;
        int pid0 = -1;
        int pid1 = -1;
        double xmin = 0.0;
        double xmax = 0.0;
    };

    auto buildSegmentsForJ =
        [&](const Geometry& geom,
            const auto& predictedDisp,
            int j,
            int nSeg,
            std::vector<ContactSeg>& out) {
            out.clear();
            out.reserve(std::max(0, nSeg - 1));

            for (int i = 1; i < nSeg; ++i) {
                const int pid0 = geom.surfp[i][j];
                const int pid1 = geom.surfp[i + 1][j];

                if (pid0 < 0 || pid1 < 0) continue;

                const double x0 = predictedDisp[pid0].ux;
                const double x1 = predictedDisp[pid1].ux;

                const double xmin = std::min(x0, x1);
                const double xmax = std::max(x0, x1);

                if (xmax - xmin <= EPS_X) continue;

                out.push_back(ContactSeg{i, pid0, pid1, xmin, xmax});
            }

            std::sort(out.begin(), out.end(),
                [](const ContactSeg& a, const ContactSeg& b) {
                    return a.xmin < b.xmin;
                });
        };

    auto processPair =
        [&](int j,
            const ContactSeg& L,
            const ContactSeg& R) {
            const int iL = L.i;
            const int iR = R.i;

            const int pidL0 = L.pid0;
            const int pidL1 = L.pid1;
            const int pidR0 = R.pid0;
            const int pidR1 = R.pid1;

            const double Lxmin = L.xmin;
            const double Lxmax = L.xmax;
            const double Rxmin = R.xmin;
            const double Rxmax = R.xmax;

            const double xa = std::max(Lxmin, Rxmin);
            const double xb = std::min(Lxmax, Rxmax);

            if (xb <= xa + EPS_X) return;
            double sLa, sRa, sLb, sRb;
            Pos pLa, pRa, pLb, pRb;

            const double gapA =
                evalPredAtX(pidL0, pidL1, pidR0, pidR1,
                            xa, sLa, sRa, pLa, pRa);

            const double gapB =
                evalPredAtX(pidL0, pidL1, pidR0, pidR1,
                            xb, sLb, sRb, pLb, pRb);

            const double pairZMm = 0.25 * (
                pLa.z + pLb.z + pRa.z + pRb.z);
            double xc0, xc1;
            double xContact;
            double penMean_mm;
            double penMax_mm;

            const bool hasPenetration =
                computePenetrationCentroid(
                    xa, xb,
                    gapA, gapB,
                    xc0, xc1,
                    xContact,
                    penMean_mm,
                    penMax_mm
                );

            if (!hasPenetration) return;

            double sL, sR;
            Pos pL, pR;

            evalPredAtX(pidL0, pidL1, pidR0, pidR1,
                        xContact, sL, sR, pL, pR);

            double nx = 0.0;
            double ny = 0.0;
            computeXYNormal(pidL0, pidL1, pidR0, pidR1, pL, pR, nx, ny);

            // Penetration and its rate use exactly the same normal as the
            // reaction force.
            const double normalGap = (pR.x - pL.x) * nx + (pR.y - pL.y) * ny;
            const double pen = std::max(0.0, -normalGap) * 1.0e-3;
            if (pen <= 0.0 || !std::isfinite(pen)) return;

            const Pos vL = interpVelocity(stateL.vel, pidL0, pidL1, sL);
            const Pos vR = interpVelocity(stateR.vel, pidR0, pidR1, sR);
            // vel is stored in mm/s; only a closing normal velocity adds
            // one-sided damping pressure.
            const double normalSeparationSpeed =
                (vR.x - vL.x) * nx + (vR.y - vL.y) * ny;
            const double penDot = std::max(0.0, -normalSeparationSpeed) * 1.0e-3;

            const double non_linear_term =
                1.0 + sp.kc2 * contactOmega2 * pen * pen;

            const double f_spring =
                sp.kc1 * contactOmega2 * pen * non_linear_term;

            const double f_damp =
                sp.kc3 * penDot;

            double contact_pressure = f_spring + f_damp;

            if (contact_pressure < 0.0) contact_pressure = 0.0;
            if (!std::isfinite(contact_pressure)) return;
            if (sp.kc1 == 0.0 && sp.kc2 == 0.0 && sp.kc3 == 0.0
                && std::abs(contact_pressure) > 1.0e-14) {
                std::cerr << "[CONTACT ERROR] zero coefficients but nonzero pressure: "
                          << contact_pressure << " Pa at step=" << step
                          << " iter=" << contactIter << "\n";
            }

            const double overlapLen = xc1 - xc0;

            const double lenXL = Lxmax - Lxmin;
            const double lenXR = Rxmax - Rxmin;

            const double fracL = clamp01(overlapLen / std::max(lenXL, EPS_X));
            const double fracR = clamp01(overlapLen / std::max(lenXR, EPS_X));

            const double areaL =
                0.5 * (geomL.sarea[iL][j] + geomL.sarea[iL + 1][j])
                * fracL * 1.0e-6;

            const double areaR =
                0.5 * (geomR.sarea[iR][j] + geomR.sarea[iR + 1][j])
                * fracR * 1.0e-6;

            const double area = 0.5 * (areaL + areaR);

            if (area <= 0.0 || !std::isfinite(area)) return;

            const double F = contact_pressure * area;
            const bool excluded = contactIsExcluded(j, pairZMm);

            const double FxL = -F * nx;
            const double FyL = -F * ny;

            const double FxR =  F * nx;
            const double FyR =  F * ny;

            if (pen > maxPenetration) {
                maxPenetration = pen;
            }

            if (!excluded) {
                addSegmentForceXY(
                    fdisXL, fdisYL,
                    iL, iL + 1, j, sL,
                    FxL, FyL
                );

                addSegmentForceXY(
                    fdisXR, fdisYR,
                    iR, iR + 1, j, sR,
                    FxR, FyR
                );
            }

            contactFlag = true;

#ifdef PROFILE_CALCDIS
            call_contacts++;
#endif
        };

    std::vector<ContactSeg> segL;
    std::vector<ContactSeg> segR;
    if (!diagnosticOnly
        && (contactIter == 1 || static_cast<int>(contactPairCache.size()) != maxJ)) {
        contactPairCache.assign(maxJ, {});
    }

    auto rangesAreOrdered = [&](const std::vector<ContactSeg>& segments) {
        for (std::size_t k = 1; k < segments.size(); ++k) {
            if (segments[k].xmin < segments[k - 1].xmax - EPS_X) return false;
        }
        return true;
    };

    for (int j = 1; j < maxJ - 1; ++j) {

#ifdef PROFILE_CALCDIS
        auto tb0 = now();
#endif

        buildSegmentsForJ(geomL, stateL.predictedDisp, j, nSegL, segL);
        buildSegmentsForJ(geomR, stateR.predictedDisp, j, nSegR, segR);

#ifdef PROFILE_CALCDIS
        {
            auto tb1 = now();
            prof_build_ms += elapsed_ms(tb0, tb1);

            call_segL += static_cast<long long>(segL.size());
            call_segR += static_cast<long long>(segR.size());
        }
#endif

        std::vector<int> indexL(nSegL + 1, -1);
        std::vector<int> indexR(nSegR + 1, -1);
        for (int k = 0; k < static_cast<int>(segL.size()); ++k) indexL[segL[k].i] = k;
        for (int k = 0; k < static_cast<int>(segR.size()); ++k) indexR[segR[k].i] = k;

        auto runCandidate = [&](const ContactSeg& L, const ContactSeg& R) {
            if (std::min(L.xmax, R.xmax) <= std::max(L.xmin, R.xmin) + EPS_X) return;
            processPair(j, L, R);
#ifdef PROFILE_CALCDIS
            ++call_process_pairs;
#endif
        };

        // Later contact iterations use the candidate topology found in the
        // first iteration.  A new time step always performs a full search.
        if (contactIter > 1 && !contactPairCache[j].empty()) {
            for (const auto& [iL, iR] : contactPairCache[j]) {
                if (iL < 0 || iL >= static_cast<int>(indexL.size()) ||
                    iR < 0 || iR >= static_cast<int>(indexR.size())) continue;
                const int l = indexL[iL], r = indexR[iR];
                if (l < 0 || r < 0) continue;
                runCandidate(segL[l], segR[r]);
            }
            continue;
        }

        if (rangesAreOrdered(segL) && rangesAreOrdered(segR)) {
            std::size_t l = 0, r = 0;
            while (l < segL.size() && r < segR.size()) {
                const auto& L = segL[l];
                const auto& R = segR[r];
                if (std::min(L.xmax, R.xmax) > std::max(L.xmin, R.xmin) + EPS_X) {
                    runCandidate(L, R);
                    if (!diagnosticOnly && contactIter == 1)
                        contactPairCache[j].emplace_back(L.i, R.i);
                }
                if (L.xmax < R.xmax - EPS_X) ++l;
                else if (R.xmax < L.xmax - EPS_X) ++r;
                else { ++l; ++r; }
            }
        } else {
            // Non-monotone x geometry is rare; retain a complete, safe
            // fallback rather than silently missing contact.
            for (const auto& L : segL) for (const auto& R : segR) {
                if (std::min(L.xmax, R.xmax) <= std::max(L.xmin, R.xmin) + EPS_X) continue;
                runCandidate(L, R);
                if (!diagnosticOnly && contactIter == 1)
                    contactPairCache[j].emplace_back(L.i, R.i);
            }
        }
    }

    if (diagnosticOnly) {
        contactFlag = savedContactFlag;
        max_force_diff = savedMaxForceDiff;
        contact_force_residual = savedForceResidual;
        contact_penetration_residual = savedPenResidual;
        max_contact_penetration = savedMaxPen;
        previous_contact_penetration_ = savedPreviousPen;
        return;
    }

#ifdef PROFILE_CALCDIS
    auto ta0 = now();
#endif

    max_force_diff = 0.0;
    double contactForceNormSquared = 0.0;

    for (int i = 0; i < static_cast<int>(fdisXL.size()); ++i) {
        const int nj = std::min({
            static_cast<int>(fdisXL[i].size()),
            static_cast<int>(fxL[i].size()),
            static_cast<int>(prevFdisXL[i].size())
        });

        for (int j = 0; j < nj; ++j) {
            const double dFx = fdisXL[i][j] - prevFdisXL[i][j];
            const double dFy = fdisYL[i][j] - prevFdisYL[i][j];

            fxL[i][j] += dFx;
            fyL[i][j] += dFy;

            max_force_diff = std::max(max_force_diff, std::abs(dFx));
            max_force_diff = std::max(max_force_diff, std::abs(dFy));
            contactForceNormSquared += fdisXL[i][j] * fdisXL[i][j]
                                    + fdisYL[i][j] * fdisYL[i][j];

            prevFdisXL[i][j] = fdisXL[i][j];
            prevFdisYL[i][j] = fdisYL[i][j];
        }
    }

    for (int i = 0; i < static_cast<int>(fdisXR.size()); ++i) {
        const int nj = std::min({
            static_cast<int>(fdisXR[i].size()),
            static_cast<int>(fxR[i].size()),
            static_cast<int>(prevFdisXR[i].size())
        });

        for (int j = 0; j < nj; ++j) {
            const double dFx = fdisXR[i][j] - prevFdisXR[i][j];
            const double dFy = fdisYR[i][j] - prevFdisYR[i][j];

            fxR[i][j] += dFx;
            fyR[i][j] += dFy;

            max_force_diff = std::max(max_force_diff, std::abs(dFx));
            max_force_diff = std::max(max_force_diff, std::abs(dFy));
            contactForceNormSquared += fdisXR[i][j] * fdisXR[i][j]
                                    + fdisYR[i][j] * fdisYR[i][j];

            prevFdisXR[i][j] = fdisXR[i][j];
            prevFdisYR[i][j] = fdisYR[i][j];
        }
    }

    max_contact_penetration = maxPenetration;
    contact_penetration_residual = std::abs(maxPenetration - previous_contact_penetration_);
    previous_contact_penetration_ = maxPenetration;
    contact_force_residual = max_force_diff /
        std::max(std::sqrt(contactForceNormSquared), 1.0e-9);

#ifdef PROFILE_CALCDIS
    {
        auto ta1 = now();
        prof_accum_ms += elapsed_ms(ta0, ta1);
    }

    {
        auto prof_t_total1 = now();
        prof_total_ms += elapsed_ms(prof_t_total0, prof_t_total1);
    }

    prof_calls++;

    prof_total_segL += call_segL;
    prof_total_segR += call_segR;
    prof_total_process_pairs += call_process_pairs;
    prof_total_contacts += call_contacts;
    prof_total_active_candidates += call_active_candidates;
    prof_max_activeR = std::max(prof_max_activeR, call_max_activeR);

    if (prof_calls % PROFILE_CALCDIS_EVERY == 0) {
        const double inv_calls = 1.0 / static_cast<double>(prof_calls);

        const double avg_segL_per_call =
            static_cast<double>(prof_total_segL) * inv_calls;

        const double avg_segR_per_call =
            static_cast<double>(prof_total_segR) * inv_calls;

        const double avg_pairs_per_call =
            static_cast<double>(prof_total_process_pairs) * inv_calls;

        const double avg_contacts_per_call =
            static_cast<double>(prof_total_contacts) * inv_calls;

        const double contact_ratio =
            static_cast<double>(prof_total_contacts) /
            static_cast<double>(std::max<long long>(1, prof_total_process_pairs));

        const double avg_active_candidates_per_call =
            static_cast<double>(prof_total_active_candidates) * inv_calls;

    }
#endif
}

void ForceCalculator::resetPreviousContactForce() {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(prevFdisXL.size()); ++i) std::fill(prevFdisXL[i].begin(), prevFdisXL[i].end(), 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(prevFdisYL.size()); ++i) std::fill(prevFdisYL[i].begin(), prevFdisYL[i].end(), 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(prevFdisXR.size()); ++i) std::fill(prevFdisXR[i].begin(), prevFdisXR[i].end(), 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(prevFdisYR.size()); ++i) std::fill(prevFdisYR[i].begin(), prevFdisYR[i].end(), 0.0);
    previous_contact_penetration_ = 0.0;
    contact_force_residual = 0.0;
    contact_penetration_residual = 0.0;
    max_contact_penetration = 0.0;

}
