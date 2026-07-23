#include "ForceCalculator.h"
#include "Displacement.h"
#include "FluidConstants.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <chrono>


#ifndef PROFILE_CALCDIS_EVERY
#define PROFILE_CALCDIS_EVERY 1000
#endif


auto checkNaN = [](double val, const std::string& name) {
    if (std::isnan(val) || std::isinf(val)) {
        std::cerr << "[NaN DETECTED] " << name << " = " << val << std::endl;
    }
};



ForceCalculator::ForceCalculator(const Geometry& geomL_, const Geometry& geomR_, 
                                 const ModeData& mdL_, const ModeData& mdR_, 
                                 State& stL_, State& stR_, const SimulationParams& sp_)
    : geomL(geomL_), geomR(geomR_), 
      modeDataL(mdL_), modeDataR(mdR_), 
      stateL(stL_), stateR(stR_), sp(sp_),
      fxL(surfaceLoad.fxL), fyL(surfaceLoad.fyL), fzL(surfaceLoad.fzL),
      fxR(surfaceLoad.fxR), fyR(surfaceLoad.fyR), fzR(surfaceLoad.fzR),
      harea(sections.area) {}

void ForceCalculator::setOutputDirectory(const std::filesystem::path& directory) {
    outputDirectory_ = directory;
    std::filesystem::create_directories(outputDirectory_);
}


void ForceCalculator::initialize() {
    // --- 左右それぞれのサイズを取得 ---
    int nsurfl_L = geomL.nsurfl;
    int nsurfz_L = geomL.nsurfz;
    int nModes_L = modeDataL.nModes;

    int nsurfl_R = geomR.nsurfl;
    int nsurfz_R = geomR.nsurfz;
    int nModes_R = modeDataR.nModes;

    // --- 左声帯 (L) の変数初期化 ---
    surfaceLoad.resize(geomL, geomR);
    fiL.assign(nModes_L, 0.0);

    fluidPressureL.assign(geomL.nsurfl, std::vector<double>(geomL.nsurfz, std::numeric_limits<double>::quiet_NaN()));
    fluidPressureR.assign(geomR.nsurfl, std::vector<double>(geomR.nsurfz, std::numeric_limits<double>::quiet_NaN()));
    fluidPlaneLowerL.assign(geomL.nsurfl, std::vector<int>(geomL.nsurfz, -1));
    fluidPlaneLowerR.assign(geomR.nsurfl, std::vector<int>(geomR.nsurfz, -1));
    fluidPlaneUpperWeightL.assign(geomL.nsurfl, std::vector<double>(geomL.nsurfz, std::numeric_limits<double>::quiet_NaN()));
    fluidPlaneUpperWeightR.assign(geomR.nsurfl, std::vector<double>(geomR.nsurfz, std::numeric_limits<double>::quiet_NaN()));

    // --- 右声帯 (R) の変数初期化 ---
    fiR.assign(nModes_R, 0.0);

    // --- 内部バッファ (流路断面数は左右で共通と仮定) ---
    nxsup = geomL.nxsup;

    fdisXL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));
    fdisYL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));

    fdisXR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));
    fdisYR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));

    prevFdisXL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));
    prevFdisYL.assign(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));

    prevFdisXR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));
    prevFdisYR.assign(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));

    // --- 共有変数の初期化 ---
    constexpr int flowSectionCount = 50;
    sections.resize(geomL, geomR, flowSectionCount);
    psurf.assign(flowSectionCount, 0.0);
    Ug.assign(stateL.nSteps, 0.0);       // ステップ数はL/R共通
    minHarea.assign(stateL.nSteps, 0.0);

    currentUg = 0.0;    
    flowModel.initialize(sp, geomL, stateL.nSteps);

    rho = sp.rho ;
    mu  = sp.mu ;
    double zmin = std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (const auto& point : geomL.points) {
        zmin = std::min(zmin, point.z);
        zmax = std::max(zmax, point.z);
    }
    lg = zmax - zmin;

    contactForceL_ij.resize(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));
    contactForceR_ij.resize(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));


    //[DEBUG]
    debugForceFile.open(outputDirectory_ / "debug_force.txt", std::ios::trunc);
    if (debugForceFile) {
        debugForceFile << "=== Debug Force Log Initialized ===\n";
    }

    contactDebugFile.open(outputDirectory_ / "contact_debug.csv", std::ios::trunc);
    if (contactDebugFile) {
        contactDebugFile
            << "step,iter,contacts,attractive_contacts,nonfinite_contacts,"
            << "max_pen_m,max_abs_pen_dot,max_contact_pressure,max_force,"
            << "min_sep_dot_norm,worst_iL,worst_iR,worst_j,worst_nx,worst_ny,"
            << "worst_FxL,worst_FyL,worst_gapC,worst_gapPrev\n"
            << "sumContactL,sumContactR,"
            << "maxContactL,maxContactIL,maxContactJL,"
            << "maxContactR,maxContactIR,maxContactJR,";            
    }

    contactMonitorFile.open(outputDirectory_ / "contact_monitor.csv", std::ios::trunc);
    if (contactMonitorFile) {
        contactMonitorFile
            << "step,time,contact_iter,monitor_i,monitor_j,"
            << "gap_y_mm,candidate_pairs,penetrating_pairs,"
            << "contact_force_L_N,contact_force_R_N,"
            << "fdisXL_N,fdisYL_N,fdisXR_N,fdisYR_N,"
            << "contact_flag,max_force_diff,force_residual,penetration_residual_m\n";
    }

    contactIterationFile.open(outputDirectory_ / "contact_iteration_debug.csv", std::ios::trunc);
    if (contactIterationFile) {
        contactIterationFile
            << "step,time,contact_iter,contact_pairs,"
            << "max_pen_m,max_pen_iL,max_pen_iR,max_pen_j,"
            << "x_contact_mm,gap_a_mm,gap_b_mm,gap_c_mm,"
            << "closing_speed_mps,contact_pressure_pa,contact_force_N,"
            << "nx,ny,normal_separation_alignment,"
            << "contact_flag,max_force_diff,force_residual,penetration_residual_m\n";
    }

    contactSearchFile.open(outputDirectory_ / "contact_search_debug.csv", std::ios::trunc);
    if (contactSearchFile) {
        contactSearchFile
            << "step,time,contact_iter,candidate_pairs,penetrating_pairs,"
            << "cached_candidates,two_pointer_candidates,fallback_candidates\n";
    }

    flowDebugFile.open(outputDirectory_ / "debug_flow_detail.csv", std::ios::trunc);
    if (flowDebugFile) {
        flowDebugFile
            << "step,time,"
            << "minA_mm2,minA_m2,"
            << "previousUg,currentUg,dUg,"
            << "currentPg,"
            << "PuLast,PuGlot,Pd0,Pd9,"
            << "nsep,"
            << "maxAbsPsurf,idxMaxAbsPsurf,"
            << "psurf0,psurfMinA,psurfLast,"
            << "pressure_recovery_clamps,first_clamp_index,first_clamp_delta_pa,"
            << "hasNonFinite\n";
    }
    
    std::cout << "[ForceCalculator] initialized (Asymmetric): "
              << "nModesL=" << nModes_L << ", nModesR=" << nModes_R
              << ", nxsup=" << nxsup 
              << ", L_sub=" << sp.L_sub
              << ", L_vt=" << sp.L_vt
              << std::endl;
    std::cout << "[ChannelSections] " << flowSectionCount
              << " fixed flow planes; retained initial constriction at i="
              << sections.referenceConstrictionI
              << " (sigma=" << sections.referenceConstrictionSigma << ")"
              << std::endl;
}

void ForceCalculator::setContactMonitor(int surfaceI, int surfaceJ) {
    contactMonitorI = surfaceI;
    contactMonitorJ = surfaceJ;
    std::cout << "[Contact] monitor surface point set to (i,j)= ("
              << contactMonitorI << ", " << contactMonitorJ << ")\n";
}

void ForceCalculator::applyFluidLoads(double t, int n) {
    // 左右で共通の分割数を使用
    int nsurfl = geomL.nsurfl;
    int nsurfz = geomL.nsurfz;
    const int nSurfaceSections = geomL.nxsup;
    const int nFlowSections = static_cast<int>(harea.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fxL.size()); ++i) {
        std::fill(fxL[i].begin(), fxL[i].end(), 0.0);
        std::fill(fyL[i].begin(), fyL[i].end(), 0.0);
        std::fill(fzL[i].begin(), fzL[i].end(), 0.0);
        std::fill(fluidPressureL[i].begin(), fluidPressureL[i].end(), std::numeric_limits<double>::quiet_NaN());
        std::fill(fluidPlaneLowerL[i].begin(), fluidPlaneLowerL[i].end(), -1);
        std::fill(fluidPlaneUpperWeightL[i].begin(), fluidPlaneUpperWeightL[i].end(), std::numeric_limits<double>::quiet_NaN());
    }

    // 右側を右側サイズでクリア
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fxR.size()); ++i) {
        std::fill(fxR[i].begin(), fxR[i].end(), 0.0);
        std::fill(fyR[i].begin(), fyR[i].end(), 0.0);
        std::fill(fzR[i].begin(), fzR[i].end(), 0.0);
        std::fill(fluidPressureR[i].begin(), fluidPressureR[i].end(), std::numeric_limits<double>::quiet_NaN());
        std::fill(fluidPlaneLowerR[i].begin(), fluidPlaneLowerR[i].end(), -1);
        std::fill(fluidPlaneUpperWeightR[i].begin(), fluidPlaneUpperWeightR[i].end(), std::numeric_limits<double>::quiet_NaN());
    }

    // ここで旧 fdisL/fdisR はクリアしない
    // 新しい fdisXL/fdisYL/fdisXR/fdisYR は calcDis() 側で管理する

    if (sp.iforce == 1) {
        // ==== sin波加振 ====
        minHarea[n] = *std::min_element(harea.begin(), harea.end());

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
        double minA = findMinHarea();
        minHarea[n] = minA;

        previousUg = (n > 0 && n-1 < (int)Ug.size()) ? Ug[n-1] : 0.0;

        // 安全な面積を使って流体計算を進める
        flowModel.advance(t, sp.dt, minA * 1e-6, previousUg, sections);
        currentUg = flowModel.flowRate();
        currentPg = flowModel.glottalPressure();
        
        // 剥離点の特定
        int nsep = findNsep(minA);

        Ug[n] = currentUg;

        //[DEBUG]
        auto maxAbsWithIndex = [](const std::vector<double>& v) {
            double maxAbs = 0.0;
            int imax = -1;
            bool bad = false;

            for (int i = 0; i < static_cast<int>(v.size()); ++i) {
                double x = v[i];
                if (!std::isfinite(x)) {
                    return std::tuple<double,int,bool>(
                        std::numeric_limits<double>::infinity(), i, true
                    );
                }
                if (std::abs(x) > maxAbs) {
                    maxAbs = std::abs(x);
                    imax = i;
                }
            }
            return std::tuple<double,int,bool>(maxAbs, imax, bad);
        };

        // Construct the pressure only up to the physical separation plane.
        // Downstream loading is mapped with a smooth physical-x blend below.
        const int sepIndex = std::clamp(nsep - 1, 0, nFlowSections - 1);
        const auto& tractPressure = flowModel.downstreamPressure();
        const double downstreamPressure = tractPressure.empty() ? 0.0 : tractPressure.front();
        int pressureRecoveryClampCount = 0;
        int firstPressureRecoveryClampIndex = -1;
        double firstPressureRecoveryDelta = 0.0;

        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = currentPg;

        if(debugForceFile){
            debugForceFile << "Step: " << std::setw(4) << n 
                           << " | minA: " << std::scientific << std::setprecision(12) << minA 
                           << " | subPg: " << std::scientific << std::setprecision(12) << currentPg  << "\n";
        }

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
                    ++pressureRecoveryClampCount;
                    if (firstPressureRecoveryClampIndex < 0) {
                        firstPressureRecoveryClampIndex = i;
                        firstPressureRecoveryDelta = psurf[i] - psurf[i - 1];
                    }
                    psurf[i] = psurf[i - 1];
                    for (int k = i + 1; k <= sepIndex; ++k) psurf[k] = psurf[i];
                    break;
                }

                //[DEBUG]
                if (n % 10 == 0 || t > 0.12) {
                    if (!std::isfinite(dx) ||
                        !std::isfinite(g3) ||
                        !std::isfinite(h_prev) ||
                        !std::isfinite(h_curr) ||
                        !std::isfinite(Ugm) ||
                        !std::isfinite(bernoulli_term) ||
                        !std::isfinite(viscous_term) ||
                        !std::isfinite(psurf[i]) ||
                        std::abs(psurf[i]) > 1.0e8) {

                        if (flowDebugFile) {
                            flowDebugFile << std::scientific << std::setprecision(12)
                                        << "#BAD_PSURF,"
                                        << "step=" << n
                                        << ",time=" << t
                                        << ",i=" << i
                                        << ",dx=" << dx
                                        << ",G3=" << g3
                                        << ",h_prev=" << h_prev
                                        << ",h_curr=" << h_curr
                                        << ",Ugm=" << Ugm
                                        << ",bernoulli=" << bernoulli_term
                                        << ",viscous=" << viscous_term
                                        << ",psurf_prev=" << psurf[i-1]
                                        << ",psurf=" << psurf[i]
                                        << "\n";
                        }
                    }
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
        struct PressureSample {
            double pressure;
            int lowerPlane;
            double upperWeight;
        };
        auto pressureAtSurfaceX = [&](double xNode) {
            if (!std::isfinite(xNode) || !std::isfinite(xSep)) {
                return PressureSample{0.0, -1, std::numeric_limits<double>::quiet_NaN()};
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
                return PressureSample{
                    (1.0 - w) * psurf[bestSegment] + w * psurf[bestSegment + 1],
                    bestSegment, w
                };
            }

                if (xNode < sections.x.front()) {
                    return PressureSample{psurf.front(), 0, 0.0};
                }
                return PressureSample{separationPressure, sepIndex, 0.0};
            }

            // A two-plane smoothstep transition suppresses a non-physical
            // force cliff after separation while tending to the tract value.
            if (std::isfinite(xBlendEnd) && xNode < xBlendEnd
                && xBlendEnd > xSep + 1.0e-12) {
                const double s = std::clamp((xNode - xSep) / (xBlendEnd - xSep), 0.0, 1.0);
                const double smooth = s * s * (3.0 - 2.0 * s);
                return PressureSample{
                    (1.0 - smooth) * separationPressure + smooth * downstreamPressure,
                    sepIndex, smooth
                };
            }
            return PressureSample{downstreamPressure, sepIndex, 1.0};
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
                    const PressureSample sample = pid_L >= 0
                        ? pressureAtSurfaceX(stateL.disp[pid_L].ux)
                        : PressureSample{0.0, -1, std::numeric_limits<double>::quiet_NaN()};
                    const double pressure = sample.pressure;
                    fluidPressureL[i][j] = pressure;
                    fluidPlaneLowerL[i][j] = sample.lowerPlane;
                    fluidPlaneUpperWeightL[i][j] = sample.upperWeight;
                    applyPressureForce(geomL, stateL, true, i, j, pressure,
                                       fxL[i][j], fyL[i][j], fzL[i][j]);

                    // if (i == 38 && j == 10) { // 情報過多を防ぐため、jは中央の代表点のみを出力
                    //     std::ofstream debugFile("../output/debug_force.txt", std::ios::app);
                    //     if(debugFile){
                    //         debugFile << std::scientific << std::setprecision(8);
                    //         debugFile << "=== Debug Step: " << n << " | i=" << i << ", j=" << j << " ===\n";
                    //         debugFile << "[Common] psurf : " << psurf[i] << "\n";
                            
                    //         debugFile << "[Left ] dsL: " << dsL << " | dzL: " << dzL 
                    //                 << " | degL0(x-y): " << degreeL[0][i][j] 
                    //                 << " | degL1(y-z): " << degreeL[1][i][j] 
                    //                 << " => fyL: " << fyL[i][j] << "\n";
                    //                 debugFile.close();}
                    // }
                }

                int pid_ip1_R = geomR.surfp[i+1][j];
                int pid_im1_R = geomR.surfp[i-1][j];
                int pid_jp1_R = geomR.surfp[i][j+1];
                int pid_jm1_R = geomR.surfp[i][j-1];

                if (pid_ip1_R >= 0 && pid_im1_R >= 0 && pid_jp1_R >= 0 && pid_jm1_R >= 0) {
                    const int pid_R = geomR.surfp[i][j];
                    const PressureSample sample = pid_R >= 0
                        ? pressureAtSurfaceX(stateR.disp[pid_R].ux)
                        : PressureSample{0.0, -1, std::numeric_limits<double>::quiet_NaN()};
                    const double pressure = sample.pressure;
                    fluidPressureR[i][j] = pressure;
                    fluidPlaneLowerR[i][j] = sample.lowerPlane;
                    fluidPlaneUpperWeightR[i][j] = sample.upperWeight;
                    applyPressureForce(geomR, stateR, false, i, j, pressure,
                                       fxR[i][j], fyR[i][j], fzR[i][j]);

                    // if (i == 38 && j == 10) { 
                    //     std::ofstream debugFile("../output/debug_force.txt", std::ios::app);
                    //     if(debugFile){
                    //         debugFile << std::scientific << std::setprecision(8);
                                    
                    //         debugFile << "[Right] dsR: " << dsR << " | dzR: " << dzR 
                    //                 << " | degR0(x-y): " << degreeR[0][i][j] 
                    //                 << " | degR1(y-z): " << degreeR[1][i][j] 
                    //                 << " => fyR: " << fyR[i][j] << "\n";
                            
                    //         debugFile << "[Diff ] abs(fyL) - abs(fyR) = " 
                    //                 << std::abs(fyL[i][j]) - std::abs(fyR[i][j]) << "\n";
                    //         debugFile << "=======================================\n";
                            
                    //         debugFile.close();}
                    // }
                }
                
                
            }
        }

        // Log only after psurf and the surface mapping have both been built.
        // Every value in one CSV row now refers to this same fluid update.
        if (flowDebugFile && (n % 10 == 0 || t > 0.12)) {
            auto [maxPsurfAbs, idxMaxPsurf, badPsurf] = maxAbsWithIndex(psurf);
            const auto& modelPu = flowModel.upstreamPressure();
            const auto& modelPd = flowModel.downstreamPressure();
            const bool hasBad = !std::isfinite(minA) || !std::isfinite(previousUg)
                || !std::isfinite(currentUg) || !std::isfinite(currentPg) || badPsurf;
            const int Nsecg = sp.N_sub;
            flowDebugFile << std::scientific << std::setprecision(12)
                << n << "," << t << "," << minA << "," << minA * 1.0e-6 << ","
                << previousUg << "," << currentUg << "," << currentUg - previousUg << ","
                << currentPg << ","
                << (Nsecg >= 0 && Nsecg < static_cast<int>(modelPu.size()) ? modelPu[Nsecg] : 0.0) << ","
                << (Nsecg + 1 < static_cast<int>(modelPu.size()) ? modelPu[Nsecg + 1] : 0.0) << ","
                << (modelPd.empty() ? 0.0 : modelPd.front()) << ","
                << (modelPd.empty() ? 0.0 : modelPd.back()) << ","
                << nsep << "," << maxPsurfAbs << "," << idxMaxPsurf << ","
                << (psurf.empty() ? 0.0 : psurf.front()) << ","
                << (idxMaxPsurf >= 0 ? psurf[idxMaxPsurf] : 0.0) << ","
                << (psurf.empty() ? 0.0 : psurf.back()) << ","
                << pressureRecoveryClampCount << "," << firstPressureRecoveryClampIndex << ","
                << firstPressureRecoveryDelta << "," << hasBad << "\n";
        }

    }
    

}

void ForceCalculator::projectLoadsToModes() {
    modalProjector.project(geomL, modeDataL, geomR, modeDataR,
                           surfaceLoad, fiL, fiR);
}

ForceCalculator::FluidVisualizationFields ForceCalculator::fluidVisualizationFields() const {
    const int nLeft = geomL.nPoints;
    const int nRight = geomR.nPoints;
    const int nTotal = nLeft + nRight;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    FluidVisualizationFields out;
    out.forceN.assign(nTotal, {0.0, 0.0, 0.0});
    out.pressurePa.assign(nTotal, nan);
    out.forceMagnitudeN.assign(nTotal, nan);
    out.planeLowerIndex.assign(nTotal, -1.0);
    out.planeUpperWeight.assign(nTotal, nan);
    out.planeCoordinate.assign(nTotal, nan);
    out.active.assign(nTotal, 0.0);

    auto transfer = [&](const Geometry& geom, int pointOffset,
                        const std::vector<std::vector<double>>& fx,
                        const std::vector<std::vector<double>>& fy,
                        const std::vector<std::vector<double>>& fz,
                        const std::vector<std::vector<double>>& pressure,
                        const std::vector<std::vector<int>>& lowerPlane,
                        const std::vector<std::vector<double>>& upperWeight) {
        const int ni = std::min(geom.nxsup, static_cast<int>(pressure.size()));
        for (int i = 0; i < ni; ++i) {
            const int nj = std::min(geom.nsurfz, static_cast<int>(pressure[i].size()));
            for (int j = 0; j < nj; ++j) {
                const int pid = geom.surfp[i][j];
                if (pid < 0 || pid >= geom.nPoints) continue;
                const int outputId = pointOffset + pid;
                const double p = pressure[i][j];
                const int k = lowerPlane[i][j];
                const double w = upperWeight[i][j];
                if (std::isfinite(p)) out.pressurePa[outputId] = p;
                if (k >= 0) {
                    out.planeLowerIndex[outputId] = static_cast<double>(k);
                    out.planeUpperWeight[outputId] = w;
                    out.planeCoordinate[outputId] = static_cast<double>(k) + w;
                    out.active[outputId] = 1.0;
                }
                if (i < static_cast<int>(fx.size()) && j < static_cast<int>(fx[i].size())) {
                    const std::array<double, 3> f = {fx[i][j], fy[i][j], fz[i][j]};
                    out.forceN[outputId] = f;
                    out.forceMagnitudeN[outputId] = std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
                }
            }
        }
    };

    transfer(geomL, 0, fxL, fyL, fzL,
             fluidPressureL, fluidPlaneLowerL, fluidPlaneUpperWeightL);
    transfer(geomR, nLeft, fxR, fyR, fzR,
             fluidPressureR, fluidPlaneLowerR, fluidPlaneUpperWeightR);
    return out;
}
// void ForceCalculator::contactForce() {

//     contactFlag = false;
//     double omg1 = 2.0 * M_PI * modeData.frequencies[0]; // 1次固有振動数
//     double omg2 = omg1 * omg1;

//     for (int i = 0; i < geom.nxsup; ++i) {           // nxsup は計算範囲
//         for (int j = 1; j < geom.nsurfz - 1; ++j) {  // 2..nsurfz-1 (0-index)
//             int node = geom.surfp[i][j];

//             double y     = state.disp[node].uy;        // 現在変位
//             double ydot  = state.vel[node].uy;         // 現在速度
//             double ymid  = geom.ymid[j];
//             double yhat  = y + sp.dt * ydot;           // 予測位置


//             // 接触状態を判定
//             bool contact_now    = (y > ymid);          // 現時点で接触
//             bool contact_future = (yhat > ymid);       // 次ステップで接触

//             if (!contact_now) {
//                 continue;
//             }

//             double pen = (ymid - y) * 1e-3; 

//             double f_contact = sp.kc1 * omg2 * pen * (1.0 + sp.kc2 * omg2 * pen * pen);

//             double f_damp =  sp.kc3 * pen * ydot;

//             double f_total = (f_contact + f_damp) * geom.sarea[i][j] * 1e-6;

//             if (f_total > 0.0) { f_total = 0.0; }           

//             fy[i][j] += f_total;
//             contactFlag = true;
//         }
//     } 

// }

void ForceCalculator::updateChannelSections() {
    sectionBuilder.build(geomL, stateL, geomR, stateR, sections);
}

void ForceCalculator::applyContactLoads(int step, int contactIter) {

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
    lineStartL.clear();
    lineEndL.clear();

    const double time = step * sp.dt;

    const bool debugContact =
        (time >= 0.120 && time <= 0.135);
    int debugContactCount = 0;
    int debugAttractiveContactCount = 0;
    int debugNonfiniteContactCount = 0;
    double debugMaxPen = 0.0;
    double debugMaxAbsPenDot = 0.0;
    double debugMaxContactPressure = 0.0;
    double debugMaxForce = 0.0;
    double debugMinSepDotNorm = std::numeric_limits<double>::infinity();
    int debugWorstIL = -1;
    int debugWorstIR = -1;
    int debugWorstJ = -1;
    double debugWorstNx = 0.0;
    double debugWorstNy = 0.0;
    double debugWorstFxL = 0.0;
    double debugWorstFyL = 0.0;
    double debugWorstGapC = 0.0;
    double debugWorstGapPrev = 0.0;
    int monitorCandidatePairs = 0;
    int monitorPenetratingPairs = 0;
    int candidatePairCount = 0;
    int cachedCandidateCount = 0;
    int twoPointerCandidateCount = 0;
    int fallbackCandidateCount = 0;
    int contactPairCount = 0;
    double maxPenetration = 0.0;
    int maxPenIL = -1, maxPenIR = -1, maxPenJ = -1;
    double maxPenX = std::numeric_limits<double>::quiet_NaN();
    double maxPenGapA = std::numeric_limits<double>::quiet_NaN();
    double maxPenGapB = std::numeric_limits<double>::quiet_NaN();
    double maxPenGapC = std::numeric_limits<double>::quiet_NaN();
    double maxPenDot = std::numeric_limits<double>::quiet_NaN();
    double maxPenPressure = std::numeric_limits<double>::quiet_NaN();
    double maxPenForce = std::numeric_limits<double>::quiet_NaN();
    double maxPenNx = std::numeric_limits<double>::quiet_NaN();
    double maxPenNy = std::numeric_limits<double>::quiet_NaN();
    double maxPenNormalAlignment = std::numeric_limits<double>::quiet_NaN();

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

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(contactForceL_ij.size()); ++i) {
        std::fill(contactForceL_ij[i].begin(), contactForceL_ij[i].end(), 0.0);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(contactForceR_ij.size()); ++i) {
        std::fill(contactForceR_ij[i].begin(), contactForceR_ij[i].end(), 0.0);
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

    auto evalPrevAtX =
        [&](int pidL0, int pidL1,
            int pidR0, int pidR1,
            double x,
            double& sL,
            double& sR,
            Pos& pL,
            Pos& pR) {
            const double xL0 = stateL.disp[pidL0].ux;
            const double xL1 = stateL.disp[pidL1].ux;
            const double xR0 = stateR.disp[pidR0].ux;
            const double xR1 = stateR.disp[pidR1].ux;

            sL = paramAtX(x, xL0, xL1);
            sR = paramAtX(x, xR0, xR1);

            pL = interpPos(stateL.disp, pidL0, pidL1, sL);
            pR = interpPos(stateR.disp, pidR0, pidR1, sR);

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
            std::vector<std::vector<double>>& contactForce,
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

            const double mag = std::sqrt(Fx * Fx + Fy * Fy);
            contactForce[i0][j] += w0 * mag;
            contactForce[i1][j] += w1 * mag;
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

    static bool is_first_contact = true;

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
            ++candidatePairCount;

            const bool touchesMonitor =
                j == contactMonitorJ &&
                ((iL <= contactMonitorI && contactMonitorI <= iL + 1) ||
                 (iR <= contactMonitorI && contactMonitorI <= iR + 1));
            if (touchesMonitor) ++monitorCandidatePairs;

            double sLa, sRa, sLb, sRb;
            Pos pLa, pRa, pLb, pRb;

            const double gapA =
                evalPredAtX(pidL0, pidL1, pidR0, pidR1,
                            xa, sLa, sRa, pLa, pRa);

            const double gapB =
                evalPredAtX(pidL0, pidL1, pidR0, pidR1,
                            xb, sLb, sRb, pLb, pRb);

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
            if (touchesMonitor) ++monitorPenetratingPairs;

            double sL, sR;
            Pos pL, pR;

            const double gapC =
                evalPredAtX(pidL0, pidL1, pidR0, pidR1,
                            xContact, sL, sR, pL, pR);

            double nx = 0.0;
            double ny = 0.0;
            computeXYNormal(pidL0, pidL1, pidR0, pidR1, pL, pR, nx, ny);

            // Penetration and its rate use exactly the same normal as the
            // reaction force.  gapC remains only a diagnostic y-gap.
            const double normalGap = (pR.x - pL.x) * nx + (pR.y - pL.y) * ny;
            const double pen = std::max(0.0, -normalGap) * 1.0e-3;
            if (pen <= 0.0 || !std::isfinite(pen)) return;
            ++contactPairCount;

            if (is_first_contact) {
                std::cout << "\n=== [DEBUG] FIRST NORMAL CONTACT DETECTED ==="
                          << std::endl;
                is_first_contact = false;
            }

            double sLprev, sRprev;
            Pos pLprev, pRprev;
            const double gapPrev =
                evalPrevAtX(pidL0, pidL1, pidR0, pidR1,
                            xContact, sLprev, sRprev, pLprev, pRprev);

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

            if (pen > maxPenetration) {
                const double sepX = pR.x - pL.x;
                const double sepY = pR.y - pL.y;
                const double sepLen = norm2(sepX, sepY);
                maxPenetration = pen;
                maxPenIL = iL;
                maxPenIR = iR;
                maxPenJ = j;
                maxPenX = xContact;
                maxPenGapA = gapA;
                maxPenGapB = gapB;
                maxPenGapC = gapC;
                maxPenDot = penDot;
                maxPenPressure = contact_pressure;
                maxPenForce = F;
                maxPenNx = nx;
                maxPenNy = ny;
                maxPenNormalAlignment = sepLen > EPS_LEN
                    ? (nx * sepX + ny * sepY) / sepLen : 0.0;
            }

            const double FxL = -F * nx;
            const double FyL = -F * ny;

            const double FxR =  F * nx;
            const double FyR =  F * ny;

            if (debugContact) {
                debugContactCount++;
                debugMaxPen = std::max(debugMaxPen, pen);
                debugMaxAbsPenDot = std::max(debugMaxAbsPenDot, std::abs(penDot));
                debugMaxContactPressure =
                    std::max(debugMaxContactPressure, contact_pressure);
                debugMaxForce = std::max(debugMaxForce, F);

                const double sepX = pL.x - pR.x;
                const double sepY = pL.y - pR.y;
                const double sepLen = std::sqrt(sepX * sepX + sepY * sepY);
                const double sepDot = FxL * sepX + FyL * sepY;
                const double sepDotNorm =
                    sepDot / std::max(F * sepLen, EPS_LEN);

                if (!std::isfinite(sepDotNorm) || !std::isfinite(F)) {
                    debugNonfiniteContactCount++;
                } else {
                    if (sepDotNorm < 0.0) {
                        debugAttractiveContactCount++;
                    }
                    if (sepDotNorm < debugMinSepDotNorm) {
                        debugMinSepDotNorm = sepDotNorm;
                        debugWorstIL = iL;
                        debugWorstIR = iR;
                        debugWorstJ = j;
                        debugWorstNx = nx;
                        debugWorstNy = ny;
                        debugWorstFxL = FxL;
                        debugWorstFyL = FyL;
                        debugWorstGapC = gapC;
                        debugWorstGapPrev = gapPrev;
                    }
                }
            }

            addSegmentForceXY(
                fdisXL, fdisYL, contactForceL_ij,
                iL, iL + 1, j, sL,
                FxL, FyL
            );

            addSegmentForceXY(
                fdisXR, fdisYR, contactForceR_ij,
                iR, iR + 1, j, sR,
                FxR, FyR
            );

            if (j == 35) {
                lineStartL.push_back({
                    stateL.predictedDisp[pidL0].ux,
                    stateL.predictedDisp[pidL0].uy,
                    stateL.predictedDisp[pidL0].uz
                });
                lineEndL.push_back({
                    stateL.predictedDisp[pidL1].ux,
                    stateL.predictedDisp[pidL1].uy,
                    stateL.predictedDisp[pidL1].uz
                });

                lineStartL.push_back({
                    stateR.predictedDisp[pidR0].ux,
                    stateR.predictedDisp[pidR0].uy,
                    stateR.predictedDisp[pidR0].uz
                });
                lineEndL.push_back({
                    stateR.predictedDisp[pidR1].ux,
                    stateR.predictedDisp[pidR1].uy,
                    stateR.predictedDisp[pidR1].uz
                });
            }

            contactFlag = true;

#ifdef PROFILE_CALCDIS
            call_contacts++;
#endif
        };

    std::vector<ContactSeg> segL;
    std::vector<ContactSeg> segR;
    if (contactIter == 1 || static_cast<int>(contactPairCache.size()) != maxJ) {
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

        auto runCandidate = [&](const ContactSeg& L, const ContactSeg& R, int& counter) {
            if (std::min(L.xmax, R.xmax) <= std::max(L.xmin, R.xmin) + EPS_X) return;
            ++counter;
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
                runCandidate(segL[l], segR[r], cachedCandidateCount);
            }
            continue;
        }

        if (rangesAreOrdered(segL) && rangesAreOrdered(segR)) {
            std::size_t l = 0, r = 0;
            while (l < segL.size() && r < segR.size()) {
                const auto& L = segL[l];
                const auto& R = segR[r];
                if (std::min(L.xmax, R.xmax) > std::max(L.xmin, R.xmin) + EPS_X) {
                    runCandidate(L, R, twoPointerCandidateCount);
                    if (contactIter == 1) contactPairCache[j].emplace_back(L.i, R.i);
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
                runCandidate(L, R, fallbackCandidateCount);
                if (contactIter == 1) contactPairCache[j].emplace_back(L.i, R.i);
            }
        }
    }

double sumContactL = 0.0;
double sumContactR = 0.0;
double maxContactL = 0.0;
double maxContactR = 0.0;
int maxContactIL = -1, maxContactJL = -1;
int maxContactIR = -1, maxContactJR = -1;

for (int i = 0; i < static_cast<int>(contactForceL_ij.size()); ++i) {
    for (int j = 0; j < static_cast<int>(contactForceL_ij[i].size()); ++j) {
        double v = contactForceL_ij[i][j];
        if (std::isfinite(v)) {
            sumContactL += v;
            if (v > maxContactL) {
                maxContactL = v;
                maxContactIL = i;
                maxContactJL = j;
            }
        }
    }
}

for (int i = 0; i < static_cast<int>(contactForceR_ij.size()); ++i) {
    for (int j = 0; j < static_cast<int>(contactForceR_ij[i].size()); ++j) {
        double v = contactForceR_ij[i][j];
        if (std::isfinite(v)) {
            sumContactR += v;
            if (v > maxContactR) {
                maxContactR = v;
                maxContactIR = i;
                maxContactJR = j;
            }
        }
    }
}

    if (debugContact && contactDebugFile) {
        if (!std::isfinite(debugMinSepDotNorm)) {
            debugMinSepDotNorm = 0.0;
        }
        contactDebugFile << step << "," << contactIter << ","
                         << debugContactCount << ","
                         << debugAttractiveContactCount << ","
                         << debugNonfiniteContactCount << ","
                         << debugMaxPen << ","
                         << debugMaxAbsPenDot << ","
                         << debugMaxContactPressure << ","
                         << debugMaxForce << ","
                         << debugMinSepDotNorm << ","
                         << debugWorstIL << ","
                         << debugWorstIR << ","
                         << debugWorstJ << ","
                         << debugWorstNx << ","
                         << debugWorstNy << ","
                         << debugWorstFxL << ","
                         << debugWorstFyL << ","
                         << debugWorstGapC << ","
                         << debugWorstGapPrev << ","
                         << sumContactL << ","
                         << sumContactR << ","
                         << maxContactL << ","
                         << maxContactIL << ","
                         << maxContactJL << ","
                         << maxContactR << ","
                         << maxContactIR << ","
                         << maxContactJR << "\n";
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

    if (contactMonitorFile && contactMonitorI >= 0 && contactMonitorJ >= 0 &&
        contactMonitorI < geomL.nxsup && contactMonitorI < geomR.nxsup &&
        contactMonitorJ < geomL.nsurfz && contactMonitorJ < geomR.nsurfz) {
        const int pidL = geomL.surfp[contactMonitorI][contactMonitorJ];
        const int pidR = geomR.surfp[contactMonitorI][contactMonitorJ];
        const double gapY = (pidL >= 0 && pidR >= 0)
            ? stateR.predictedDisp[pidR].uy - stateL.predictedDisp[pidL].uy
            : std::numeric_limits<double>::quiet_NaN();
        contactMonitorFile << std::scientific << std::setprecision(12)
            << step << "," << time << "," << contactIter << ","
            << contactMonitorI << "," << contactMonitorJ << ","
            << gapY << ","
            << monitorCandidatePairs << "," << monitorPenetratingPairs << ","
            << contactForceL_ij[contactMonitorI][contactMonitorJ] << ","
            << contactForceR_ij[contactMonitorI][contactMonitorJ] << ","
            << fdisXL[contactMonitorI][contactMonitorJ] << ","
            << fdisYL[contactMonitorI][contactMonitorJ] << ","
            << fdisXR[contactMonitorI][contactMonitorJ] << ","
            << fdisYR[contactMonitorI][contactMonitorJ] << ","
            << contactFlag << "," << max_force_diff << ","
            << contact_force_residual << "," << contact_penetration_residual << "\n";
    }

    if (contactIterationFile) {
        contactIterationFile << std::scientific << std::setprecision(12)
            << step << "," << time << "," << contactIter << ","
            << contactPairCount << ","
            << maxPenetration << ","
            << maxPenIL << "," << maxPenIR << "," << maxPenJ << ","
            << maxPenX << "," << maxPenGapA << "," << maxPenGapB << "," << maxPenGapC << ","
            << maxPenDot << "," << maxPenPressure << "," << maxPenForce << ","
            << maxPenNx << "," << maxPenNy << "," << maxPenNormalAlignment << ","
            << contactFlag << "," << max_force_diff << ","
            << contact_force_residual << "," << contact_penetration_residual << "\n";
    }

    if (contactSearchFile) {
        contactSearchFile << step << "," << std::scientific << std::setprecision(12)
            << time << "," << contactIter << ","
            << candidatePairCount << "," << contactPairCount << ","
            << cachedCandidateCount << "," << twoPointerCandidateCount << ","
            << fallbackCandidateCount << "\n";
    }


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

        // std::cout << "\n";
        // std::cout << "============================================================\n";
        // std::cout << "[PROFILE_CALCDIS] calcDis profiling summary\n";
        // std::cout << "------------------------------------------------------------\n";
        // std::cout << "calls                         : " << prof_calls << "\n";
        // std::cout << std::fixed << std::setprecision(3);
        // std::cout << "total time [ms]                : " << prof_total_ms << "\n";
        // std::cout << "  clear  [ms]                  : " << prof_clear_ms  << "\n";
        // std::cout << "  build  [ms]                  : " << prof_build_ms  << "\n";
        // std::cout << "  active [ms]                  : " << prof_active_ms << "\n";
        // std::cout << "  pair   [ms]                  : " << prof_pair_ms   << "\n";
        // std::cout << "  accum  [ms]                  : " << prof_accum_ms  << "\n";
        // std::cout << "------------------------------------------------------------\n";
        // std::cout << "avg per call [ms]              : " << prof_total_ms * inv_calls << "\n";
        // std::cout << "  clear/call  [ms]             : " << prof_clear_ms  * inv_calls << "\n";
        // std::cout << "  build/call  [ms]             : " << prof_build_ms  * inv_calls << "\n";
        // std::cout << "  active/call [ms]             : " << prof_active_ms * inv_calls << "\n";
        // std::cout << "  pair/call   [ms]             : " << prof_pair_ms   * inv_calls << "\n";
        // std::cout << "  accum/call  [ms]             : " << prof_accum_ms  * inv_calls << "\n";
        // std::cout << "------------------------------------------------------------\n";
        // std::cout << "avg segL per call              : " << avg_segL_per_call << "\n";
        // std::cout << "avg segR per call              : " << avg_segR_per_call << "\n";
        // std::cout << "avg active candidates per call : " << avg_active_candidates_per_call << "\n";
        // std::cout << "max activeR size               : " << prof_max_activeR << "\n";
        // std::cout << "processPair total calls        : " << prof_total_process_pairs << "\n";
        // std::cout << "processPair avg per call       : " << avg_pairs_per_call << "\n";
        // std::cout << "contacts total                 : " << prof_total_contacts << "\n";
        // std::cout << "contacts avg per call          : " << avg_contacts_per_call << "\n";
        // std::cout << "contact ratio                  : " << contact_ratio << "\n";
        // std::cout << "============================================================\n";
        // std::cout << "\n";
    }
#endif
}

// 1. 接触力バッファをゼロ初期化
// 2. 接触ばねに使う代表角振動数を計算
// 3. 補間・ギャップ評価・法線計算用の小関数を用意
// 4. j 断面ごとに処理
// 7. 左右セグメントの x 範囲が重なるか判定
// 8. 重なる区間の両端で y ギャップを評価
// 9. gap < 0 のめり込み区間を求める
// 10. めり込み平均量・接触代表点を求める
// 11. 前回位置との差からめり込み速度を求める
// 12. 非線形ばね + ダンパで接触圧を計算
// 13. 接触区間長と表面面積から有効面積を計算
// 14. 接触圧 × 面積で接触力を計算
// 15. XY 平面内の接触法線を計算
// 16. 左右に逆向きの力を分配
// 17. 最後に fdis を fx/fy に加算

#if 0 // Removed legacy duplicate airway solver; FlowModel is the sole owner.
void ForceCalculator::calcFlowStep(double t, double dt, double min_area) {

    static bool printed_constants = false;
    // if (t > 0.0 && !printed_constants) {
    //     std::cout << "\n=== [DEBUG] calcFlowStep Constants at Step 1 ===" << std::scientific << std::setprecision(12) << std::endl;
    //     std::cout << "ps (Lung Press): " << sp.ps << std::endl;
    //     std::cout << "rho: " << rho << " | mu: " << mu << std::endl;
    //     std::cout << "xsup: " << geomL.xsup << " | lg: " << lg << std::endl;
    //     std::cout << "beta: " << beta << " | R2: " << R2 << std::endl;
    //     std::cout << "Cu: " << Cu << " | Lu: " << Lu << std::endl;
    //     std::cout << "La: " << La << " | Ca: " << Ca << std::endl;
    //     std::cout << "=================================================\n" << std::endl;
    //     printed_constants = true;
    // }
    
    // --- 1. 声門下 (Subglottal) の更新 ---

    int Nsecp   = sp.N_vt; 
    int Nsecg   = sp.N_sub; 
    
    double rampDuration = 0.05; // 50msかけて立ち上げる
    double rampFactor = 1.0;
    
    if (t < rampDuration) {
        // Cosine Ramp (滑らか)
        rampFactor = 0.5 * (1.0 - std::cos(M_PI * t / rampDuration));
    }

    // --- ランプ適用 ---
    double currentLungPressure = sp.ps * rampFactor;

    double ug = currentUg;
    // Pu[1]...Pu[Nsecg]
    for (int j = 0; j < Nsecg; ++j) {

        Pu[j] += (Uu[j] - Uu[j+1]);
        //Pu[j] += (dt / C_use) * (Uu[j] - Uu[j+1]);
    }
    
    // 声門直下の圧力ノード (境界)
    Pu[Nsecg] +=  (Uu[Nsecg] - previousUg);
    Pu[Nsecg+1] +=  ( previousUg - Ud[0]);


    Uu[0]  -= dt / Lui * ( (dt / Cui * Pu[0]) - currentLungPressure );

    // Fortran: Uu(2)=Uu(2)-dt/(Lui+Lu)*(dt/Cu*Pu(2)-dt/Cui*Pu(1)+R2*Uu(2))
    Uu[1] -= (dt / (Lui + Lu)) * ( dt / Cu * Pu[1] - dt / Cui * Pu[0] + R2 * Uu[1] );

    for (int j = 2; j < Nsecg + 1; ++j) {
        Uu[j] -= dt / (2.0 * Lu) * ( dt / Cu * Pu[j] - dt / Cu * Pu[j-1] + R2 * Uu[j] );
    }


    // --- 2. 声門部 (Glottal Flow) の更新 ---
    if (min_area > 1e-8) {
        double min_area_m2 = min_area;
        double lis = geomL.xsup * 1e-3; // 仮定値 (Fortran側でd1+d2に相当するか要確認)
        double lg_m = lg * 1e-3;

        double Lg1 = rho *  0.5 * lis / min_area_m2;
        double Rk1 = beta * rho / ( min_area_m2 * min_area_m2); // Bernoulli (係数調整)

        double Rv1 = 12.0 * mu * lis * lg_m * lg_m / pow(min_area_m2, 3.0);

        // 駆動圧: 声門直下(Pu[last]) - 声道入口(Pd[0])
        double Ug_old = previousUg; 
        double Ug_guess = currentUg; 

        // Newton-Raphson
        for(int k=0; k<100; ++k) { // ループ回数Fortranは100
            // F(Ug)
            double F = Rk1 * std::abs(Ug_guess) * Ug_guess
                    + Rv1 * Ug_guess
                    + (Lg1 + La + Lu) * (Ug_guess - Ug_old) / dt
                    + (dt / Ca) * Pu[Nsecg+1]
                    - (dt / Cu) * Pu[Nsecg];
            
            // F'(Ug)
            double Fd =  2.0 * Rk1 * std::abs(Ug_guess) + Rv1 + (Lg1 + La + Lu) / dt;
            
            if(std::abs(F) < 1e-9) break;
            Ug_guess -= F / Fd;
        }
        currentUg = Ug_guess;
    } else {
        currentUg = 0.0;
    }

    currentPg = dt / Cu * Pu[Nsecg];


    // --- 3. 声道 (Vocal Tract) の更新 ---
    // 圧力更新 Pd[0]...

    if(!hasVocalTract ){
        for (int i = 0; i < Nsecp; i++) {
            Pd[i] = 0.0;     // 圧力ゼロ
            Ud[i] = currentUg; // 流量はすべて Ug と同じ
        }
    }else{
        Pd[0] += (dt / Ca) * (currentUg - Ud[0]);
        for(int i=1; i<Nsecp; ++i) {
            Pd[i] += (dt / Ca) * (Ud[i-1] - Ud[i]);
        }

        // 流量更新 Ud[0]...
        for(int i=0; i<Nsecp-1; ++i) {
            Ud[i] += (dt / La) * (Pd[i] - Pd[i+1]);
        }
    }
    
    // 放射端 (Radiation)
    // Lr * dUd/dt + Rr * Ud = Pd[last] - P_atm(0)
    // 離散化: (Lr/dt + Rr) * Ud_new = Pd[last] + (Lr/dt)*Ud_old
    double Z_rad = (La+Lr)/dt + Rr;
    Ud[Nsecp-1] = (Pd[Nsecp-1] + ((La+Lr)/dt)*Ud[Nsecp-1]) / Z_rad;
}
#endif

double ForceCalculator::findMinHarea() {
    double minimum = std::numeric_limits<double>::infinity();
    for (int i = 1; i + 1 < static_cast<int>(harea.size()); ++i) {
        if (!sections.valid[i] || !std::isfinite(harea[i])) continue;
        minimum = std::min(minimum, harea[i]);
    }
    return std::isfinite(minimum) ? minimum : 0.0;
}

int ForceCalculator::findNsep(double minH) {
    const int nFlowSections = static_cast<int>(harea.size());
    for (int i = 1; i + 1 < nFlowSections; i++) {
        if (sections.valid[i] && (std::fabs(harea[i] - minH) < 1e-8 || harea[i] <= 0.0)) {
            return i + 1;
        }
    }
    return nFlowSections - 1;
}



void ForceCalculator::outputForceVectors(int step) const {
    std::ostringstream stepStr;
    stepStr << std::setw(4) << std::setfill('0') << step;

    const auto forceDir = outputDirectory_ / "force_vectors";
    std::filesystem::create_directories(forceDir);
    std::ofstream fout(forceDir / ("force_" + stepStr.str() + ".csv"));
    fout << "x,y,z,Fx,Fy,Fz\n";  // CSVヘッダー

    for (int i = 1; i < geomR.nxsup - 1; ++i) {
        for (int j = 1; j < geomR.nsurfz - 1; ++j) {
            int pid = geomR.surfp[i][j];
            if (pid < 0) continue;

            const auto &p = geomR.points[pid];
            fout << p.x << "," << p.y << "," << p.z << ","
                 << fxL[i][j] << "," << fyL[i][j] << "," << fzL[i][j] << "\n";
        }
    }

    //std::cout << "[Output] force vectors written for step " << step << std::endl;
}

void ForceCalculator::outputCorrespondenceOffsets(int step) const {
    const auto offsetDir = outputDirectory_ / "pair_offset";
    std::filesystem::create_directories(offsetDir);

    // ===== 見たい対応点を指定 =====
    // 必要に応じて固定値に変えてください
    const int i_target = std::min(geomL.nxsup, geomR.nxsup) / 2;
    const int j_target = std::min(geomL.nsurfz, geomR.nsurfz) / 2;

    const int ni = std::min(geomL.nxsup, geomR.nxsup);
    const int nj = std::min(geomL.nsurfz, geomR.nsurfz);

    if (i_target < 0 || i_target >= ni || j_target < 0 || j_target >= nj) return;

    const int pidL = geomL.surfp[i_target][j_target];
    const int pidR = geomR.surfp[i_target][j_target];

    if (pidL < 0 || pidR < 0) return;

    // ===== 初期座標 =====
    const auto& pL0 = geomL.points[pidL];
    const auto& pR0 = geomR.points[pidR];

    // ===== 現在座標（変形後）=====
    const auto& pL = stateL.disp[pidL];
    const auto& pR = stateR.disp[pidR];

    // ===== 予測座標（接触判定で使っている値も見たいとき）=====
    const auto& pLf = stateL.predictedDisp[pidL];
    const auto& pRf = stateR.predictedDisp[pidR];

    // ===== 左右それぞれの変位量（本命）=====
    const double xDispL = pL.ux - pL0.x;
    const double zDispL = pL.uz - pL0.z;
    const double xDispR = pR.ux - pR0.x;
    const double zDispR = pR.uz - pR0.z;

    // ===== 参考：予測変位 =====
    const double xPredDispL = pLf.ux - pL0.x;
    const double zPredDispL = pLf.uz - pL0.z;
    const double xPredDispR = pRf.ux - pR0.x;
    const double zPredDispR = pRf.uz - pR0.z;

    // ===== 相対ズレ =====
    const double dxPair = pR.ux - pL.ux;
    const double dzPair = pR.uz - pL.uz;

    // ===== 接触判定用ギャップ =====
    const double gapY = pR.uy - pL.uy;

    // ===== CSV出力（時系列）=====
    std::ofstream fout(offsetDir / "monitor_point.csv", std::ios::app);

    // ヘッダ（初回のみ）
    if (fout.tellp() == 0) {
        fout << "step,time,"
             << "xDispL,zDispL,xDispR,zDispR,"
             << "xPredDispL,zPredDispL,xPredDispR,zPredDispR,"
             << "dxPair,dzPair,gapY\n";
    }

    // 時刻計算（dtを使う）
    double t = step * sp.dt;

    fout << step << "," << t << ","
         << xDispL << "," << zDispL << ","
         << xDispR << "," << zDispR << ","
         << xPredDispL << "," << zPredDispL << ","
         << xPredDispR << "," << zPredDispR << ","
         << dxPair << "," << dzPair << ","
         << gapY << "\n";

    fout.close();
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
