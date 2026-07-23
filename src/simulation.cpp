
#define _USE_MATH_DEFINES  
#include "Simulation.h"
#include "wavwrite.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <filesystem>




Simulation::Simulation()
    : fCalc(geomL, geomR, mdataL, mdataR, stateL, stateR, params)
{}

void Simulation::initialize(const fs::path& parameterFile) {
    std::cout << "[Simulation] Initializing..." << std::endl;
    std::string err ="error";

    if (!params.loadFromFile(parameterFile, err)) {
        throw std::runtime_error("Failed to load parameter file '" + parameterFile.string()
                                 + "': " + err);
    }
    if (!params.validate(err)) {
        throw std::runtime_error("Invalid parameter file '" + parameterFile.string()
                                 + "': " + err);
    }
    std::cout << "[Simulation] Parameter file: " << parameterFile << "\n";
    params.print();

    // Keep the normal output path stable: the analysis scripts in tools/
    // intentionally consume output/*.dat without selecting a run directory.
    // This is a latest-result workspace and is overwritten at each execution;
    // archival copies are an explicit user action rather than an unbounded
    // accumulation of heavy simulation output.
    const fs::path absoluteParameterFile = fs::absolute(parameterFile);
    const fs::path projectRoot = absoluteParameterFile.parent_path().parent_path();
    runDir = projectRoot / "output";
    fs::create_directories(runDir);
    fs::copy_file(absoluteParameterFile, runDir / "params_used.txt",
                  fs::copy_options::overwrite_existing);
    std::ofstream manifest(runDir / "manifest.txt");
    manifest << "parameter_file = " << absoluteParameterFile << "\n"
             << "nstep = " << params.nstep << "\n"
             << "dt_s = " << params.dt << "\n"
             << "nmode = " << params.nmode << "\n"
             << "ncont = " << params.ncont << "\n"
             << "iforce = " << params.iforce << "\n"
             << "contact_reference_frequency_hz = " << params.contactReferenceFrequencyHz << "\n"
             << "flow_blend_length_mm = " << params.flowBlendLengthMm << "\n"
             << "left_mode_vtu = ../input/M5_test/M5_mode_T3_b8c3.vtu\n"
             << "right_mode_vtu = ../input/M5_test/M5_mode_T3_b8c3.vtu\n"
             << "left_frequency = ../input/M5_test/M5_freq_T3_d2_b8c3.txt\n"
             << "right_frequency = ../input/M5_test/M5_freq_T3_d2_b8c3.txt\n"
             << "flow_sections = 50\n"
             << "area_close_m2 = 1e-8\n";
    fCalc.setOutputDirectory(runDir);
    std::cout << "[Simulation] Latest-result directory: " << runDir << "\n";

    geomL.loadFromVTK("../input/M5_test/M5_mode_T3_b8c3.vtu");
    geomL.surfExtractFromNAS("../input/M5_test/M5_surface_T3_d2.nas",69,70);
    geomL.surfArea();
    geomL.print();
    geomL.jtypes[5] = 3;   // 三角形
    geomL.jtypes[9] = 4;   // 四角形
    geomL.jtypes[10] = 4;
    geomL.jtypes[13] = 6;  // 六面体

    mdataL.initialize(params.nmode, geomL);

    mdataL.loadFromVTU("../input/M5_test/M5_mode_T3_b8c3.vtu", geomL);
    mdataL.loadFreqDamping("../input/M5_test/M5_freq_T3_d2_b8c3.txt");

    mdataL.normalizeModes( params.mass, geomL);
    stateL.initialize(geomL.nPoints, params.nmode, params.nstep, geomL);

    geomR.loadFromVTK("../input/M5_test/M5_mode_T3_b8c3.vtu");
    geomR.surfExtractFromNAS("../input/M5_test/M5_surface_T3_d2.nas",69,70);
    geomR.surfArea();

    geomR.jtypes[5] = 3;   // 三角形
    geomR.jtypes[9] = 4;   // 四角形
    geomR.jtypes[10] = 4;
    geomR.jtypes[13] = 6;  // 六面体

    mdataR.initialize(params.nmode, geomR);

    mdataR.loadFromVTU("../input/M5_test/M5_mode_T3_b8c3.vtu", geomR);
    mdataR.loadFreqDamping("../input/M5_test/M5_freq_T3_d2_b8c3.txt");

    mdataR.normalizeModes( params.mass, geomR);

    double ymid = geomR.ymid[0]; 

    // ① ジオメトリのYを反転
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < geomR.nPoints; ++i) {
        geomR.points[i].y = 2.0 * ymid - geomR.points[i].y;
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for (int m = 0; m < mdataR.nModes; ++m) {
        for (int i = 0; i < geomR.nPoints; ++i) {
            mdataR.modes[m][i].uy = -mdataR.modes[m][i].uy; 
        }
    }

    stateR.initialize(geomR.nPoints, params.nmode, params.nstep, geomR);

    omegaL.resize(mdataL.nModes);
    omegaR.resize(mdataR.nModes);
    for (int i = 0; i < mdataL.nModes; ++i)
        omegaL[i] = 2.0 * M_PI * mdataL.frequencies[i];
    for (int i = 0; i < mdataR.nModes; ++i)
        omegaR[i] = 2.0 * M_PI * mdataR.frequencies[i];

    fCalc.initialize();

    std::cout << "[Simulation] Initialization complete." << std::endl;
}

void Simulation::run() {
    std::cout << "[Simulation] Running..." << std::endl;

    double time_calcArea = 0.0;
    double time_calcForce = 0.0;
    double time_f2mode = 0.0;
    double time_mode2uf = 0.0;
    double time_calcDis = 0.0;
    double time_output = 0.0;

    //DEBUG
    long long total_contact_iterations = 0;
    long long steps_no_contact_break = 0;
    long long steps_converged_break = 0;
    long long steps_reached_ncont = 0;

    long long total_mode2uf_calls = 0;
    long long total_f2mode_calls = 0;
    long long total_calcDis_calls = 0;

    int max_icont_used = 0;

    auto now = []() {
        return std::chrono::high_resolution_clock::now();
    };

    auto elapsed_ms = [](auto t0, auto t1) {
        return std::chrono::duration<double, std::milli>(t1 - t0).count();
    };


    // nSteps+1 に対応
    stateL.qf.resize(mdataL.nModes, 0.0);
    stateL.qfdot.resize(mdataL.nModes, 0.0);
    stateL.qfddot.resize(mdataL.nModes, 0.0);

    stateR.qf.resize(mdataR.nModes, 0.0);
    stateR.qfdot.resize(mdataR.nModes, 0.0);
    stateR.qfddot.resize(mdataR.nModes, 0.0);

    int num = 0;

    std::ofstream fa(runDir / "area.dat");
    std::ofstream fu(runDir / "displace.dat");
    std::ofstream fp(runDir / "pressure.dat");
    std::ofstream fpv(runDir / "pressure_vt.dat");
    std::ofstream fuv(runDir / "airflow_vt.dat");
    std::ofstream fsectionx(runDir / "section_x.dat");
    std::ofstream fgapcubed(runDir / "gap_cubed.dat");
    std::ofstream fseparation(runDir / "separation.dat");
    std::ofstream fdispXY(runDir / "displace_xy.dat");
    std::ofstream fmodal(runDir / "modal_contribution.csv");
    std::ofstream fmodalDominant(runDir / "modal_dominant.csv");
    //[DEBUG]
    std::ofstream fstepdbg(runDir / "debug_step_summary.csv");
    fstepdbg << "step,time,"
            << "minArea,maxArea,idxMinArea,"
            << "currentUg,currentPg,Pd0,Pd9,"
            << "maxAbsPsurf,"
            << "maxFiL,imaxFiL,maxFiR,imaxFiR,"
            << "maxQL,maxQdL,maxQaL,maxQR,maxQdR,maxQaR,"
            << "maxPredDispL,maxPredDispR,"
            << "icont_used,contactFlag,max_force_diff,"
            << "diverged\n";
    
    std::ofstream fmodedbg(runDir / "debug_mode_summary.csv");
    fmodedbg << "step,time,icont,stage,"
            << "maxFiL,imaxFiL,maxFiR,imaxFiR,"
            << "contactFlag,max_force_diff\n";

    //std::ofstream fdbg("../output/debug_fluid.dat");
    //std::ofstream contactIterDbg("../output/contact_iteration_debug.csv");

    fa << "# step  area[mm^2]\n";
    fu << "# time_s  uyL_mm  uyR_mm\n";
    fp << "# step  pressure[Pa]\n"; 
    fpv << "# step  outlet_pressure[Pa]\n";
    fuv << "# step  airflow[m^3/s]\n";
    fsectionx << "# step  section_x[mm]\n";
    fgapcubed << "# step  integral_g_positive_cubed[mm^4]\n";
    fseparation << "# step sep_index x_sep[mm] x_blend_end[mm] p_sep[Pa]\n";
    fdispXY << "# time[s] uxL[mm] uyL[mm] uxR[mm] uyR[mm]\n";
    fmodal << "step,time_s,side,mode_index,frequency_hz,q,qdot,"
           << "probe_ux_mm,probe_uy_mm,probe_uz_mm,"
           << "surface_rms_ux_mm,surface_rms_uy_mm,surface_rms_uz_mm,"
           << "surface_rms_norm_mm\n";
    fmodalDominant << "step,time_s,side,"
                   << "dominant_probe_uy_mode,dominant_probe_uy_mm,"
                   << "dominant_surface_mode,dominant_surface_rms_mm\n";


    const double monitorTargetX = 7.2;
    const double monitorTargetZ = 8.5;
    double minDist2 = 1e100;
    int nearestIdxL = -1;
    int nearestIdxR = -1;
    int monitorI = -1;
    int monitorJ = -1;
    const int monitorNi = std::min(geomL.nxsup, geomR.nxsup);
    const int monitorNj = std::min(geomL.nsurfz, geomR.nsurfz);
    for (int i = 0; i < monitorNi; ++i) {
        for (int j = 0; j < monitorNj; ++j) {
            int idx = geomL.surfp[i][j];
            if (idx < 0 || geomR.surfp[i][j] < 0) continue;
            double dx = geomL.points[idx].x - monitorTargetX;
            double dz = geomL.points[idx].z - monitorTargetZ;
            double dist2 = dx*dx + dz*dz ;


            if (dist2 < minDist2) {
                minDist2 = dist2;
                nearestIdxL = idx;
                nearestIdxR = geomR.surfp[i][j];
                monitorI = i;
                monitorJ = j;
            }
        }
    }
    std::cout << "Monitor surface point (i,j)=(" << monitorI << ", " << monitorJ << ")\n";
    std::cout<<"Monitor Node L idx="<<geomL.points[nearestIdxL].x<<", "<<geomL.points[nearestIdxL].y<<", "<<geomL.points[nearestIdxL].z<<"\n";
    std::cout<<"Monitor Node R idx="<<geomR.points[nearestIdxR].x<<", "<<geomR.points[nearestIdxR].y<<", "<<geomR.points[nearestIdxR].z<<"\n";
    fCalc.setContactMonitor(monitorI, monitorJ);

    struct SurfaceModeRms {
        double ux = 0.0, uy = 0.0, uz = 0.0, norm = 0.0;
    };
    auto precomputeSurfaceModeRms = [](const Geometry& geom, const ModeData& modes) {
        std::vector<SurfaceModeRms> result(modes.nModes);
        for (int m = 0; m < modes.nModes; ++m) {
            double sx = 0.0, sy = 0.0, sz = 0.0;
            int count = 0;
            for (int i = 0; i < geom.nsurfl; ++i) {
                for (int j = 0; j < geom.nsurfz; ++j) {
                    const int pid = geom.surfp[i][j];
                    if (pid < 0) continue;
                    const auto& phi = modes.modes[m][pid];
                    sx += phi.ux * phi.ux;
                    sy += phi.uy * phi.uy;
                    sz += phi.uz * phi.uz;
                    ++count;
                }
            }
            if (count > 0) {
                result[m].ux = std::sqrt(sx / count);
                result[m].uy = std::sqrt(sy / count);
                result[m].uz = std::sqrt(sz / count);
                result[m].norm = std::sqrt(result[m].ux * result[m].ux
                                          + result[m].uy * result[m].uy
                                          + result[m].uz * result[m].uz);
            }
        }
        return result;
    };
    const auto surfaceModeRmsL = precomputeSurfaceModeRms(geomL, mdataL);
    const auto surfaceModeRmsR = precomputeSurfaceModeRms(geomR, mdataR);

    auto writeModalDiagnostics = [&](int step, double time, const char* side,
                                     const ModeData& modes, const State& state,
                                     int probeId,
                                     const std::vector<SurfaceModeRms>& surfaceRms) {
        int dominantProbeMode = -1;
        int dominantSurfaceMode = -1;
        double dominantProbeMagnitude = -1.0;
        double dominantSurfaceMagnitude = -1.0;
        for (int m = 0; m < modes.nModes; ++m) {
            const double scaleMm = state.q[m] * 1.0e3;
            const auto& phi = modes.modes[m][probeId];
            const double probeUx = scaleMm * phi.ux;
            const double probeUy = scaleMm * phi.uy;
            const double probeUz = scaleMm * phi.uz;
            const double rmsUx = std::abs(scaleMm) * surfaceRms[m].ux;
            const double rmsUy = std::abs(scaleMm) * surfaceRms[m].uy;
            const double rmsUz = std::abs(scaleMm) * surfaceRms[m].uz;
            const double rmsNorm = std::abs(scaleMm) * surfaceRms[m].norm;
            fmodal << std::scientific << std::setprecision(12)
                   << step << "," << time << "," << side << "," << m << ","
                   << modes.frequencies[m] << "," << state.q[m] << "," << state.qdot[m] << ","
                   << probeUx << "," << probeUy << "," << probeUz << ","
                   << rmsUx << "," << rmsUy << "," << rmsUz << "," << rmsNorm << "\n";
            if (std::abs(probeUy) > dominantProbeMagnitude) {
                dominantProbeMagnitude = std::abs(probeUy);
                dominantProbeMode = m;
            }
            if (rmsNorm > dominantSurfaceMagnitude) {
                dominantSurfaceMagnitude = rmsNorm;
                dominantSurfaceMode = m;
            }
        }
        fmodalDominant << std::scientific << std::setprecision(12)
                       << step << "," << time << "," << side << ","
                       << dominantProbeMode << "," << dominantProbeMagnitude << ","
                       << dominantSurfaceMode << "," << dominantSurfaceMagnitude << "\n";
    };


    stateL.mode2uf(geomL, mdataL, 0); 
    stateL.uf2u(); 
    stateR.mode2uf(geomR, mdataR, 0); 
    stateR.uf2u();

    std::vector<double> soundSignal;
        soundSignal.reserve(params.nstep);

    auto maxAbsVector = [](const std::vector<double>& values) {
        double maxAbs = 0.0;
        for (double v : values) {
            if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
            maxAbs = std::max(maxAbs, std::abs(v));
        }
        return maxAbs;
    };

    auto maxAbsNodeDisp = [](const Geometry& geom, const State& state) {
        double maxAbs = 0.0;
        for (int i = 0; i < state.nPoints; ++i) {
            const double dx = state.predictedDisp[i].ux - geom.points[i].x;
            const double dy = state.predictedDisp[i].uy - geom.points[i].y;
            const double dz = state.predictedDisp[i].uz - geom.points[i].z;
            if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz)) {
                return std::numeric_limits<double>::infinity();
            }
            maxAbs = std::max(maxAbs, std::abs(dx));
            maxAbs = std::max(maxAbs, std::abs(dy));
            maxAbs = std::max(maxAbs, std::abs(dz));
        }
        return maxAbs;
    };

    auto minMaxIndex = [](const std::vector<double>& v) {
        double vmin = std::numeric_limits<double>::infinity();
        double vmax = -std::numeric_limits<double>::infinity();
        int imin = -1;
        int imax = -1;

        for (int i = 0; i < static_cast<int>(v.size()); ++i) {
            double x = v[i];
            if (!std::isfinite(x)) {
                return std::tuple<double,double,int,int>(
                    -std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    i, i
                );
            }
            if (x < vmin) { vmin = x; imin = i; }
            if (x > vmax) { vmax = x; imax = i; }
        }
        return std::tuple<double,double,int,int>(vmin, vmax, imin, imax);
    };

    auto maxAbsIndex = [](const std::vector<double>& v) {
        double maxAbs = 0.0;
        int imax = -1;

        for (int i = 0; i < static_cast<int>(v.size()); ++i) {
            double x = v[i];
            if (!std::isfinite(x)) {
                return std::pair<double,int>(
                    std::numeric_limits<double>::infinity(), i
                );
            }
            if (std::abs(x) > maxAbs) {
                maxAbs = std::abs(x);
                imax = i;
            }
        }
        return std::pair<double,int>(maxAbs, imax);
    };

    auto maxAbsPsurf = [&]() {
        double m = 0.0;
        for (double p : fCalc.psurf) {
            if (!std::isfinite(p)) {
                return std::numeric_limits<double>::infinity();
            }
            m = std::max(m, std::abs(p));
        }
        return m;
    };

    //writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 1);
    //num++;
    std::cout << "[Simulation] Output step 0 (Initial State)." << std::endl;

    for (int n = 0; n < params.nstep; n++) {
        double t = n * params.dt;

        // 面積・角度の更新 (左右の相対距離で計算)
{        auto t0 = now();
        fCalc.updateChannelSections();
        auto t1 = now();
        time_calcArea += elapsed_ms(t0, t1);}


        // 圧力の計算と、左右への力の分配
{        auto t0 = now();
        fCalc.applyFluidLoads(t, n);
        auto t1 = now();
        time_calcForce += elapsed_ms(t0, t1);}

        if (n % 5 == 0) {
            fa << std::setw(4) << n;
            fp << std::setw(4) << n;
            fsectionx << std::setw(4) << n;
            fgapcubed << std::setw(4) << n;
            for (int i = 0; i < static_cast<int>(fCalc.harea.size()); ++i) {
                fa << " " << std::setw(8) << fCalc.harea[i] << " ";
                fp << " " << std::setw(8) << fCalc.psurf[i] << " ";
                fsectionx << " " << std::setw(8) << fCalc.sectionX(i) << " ";
                fgapcubed << " " << std::setw(8) << fCalc.sectionGapCubed(i) << " ";
            }
            fa << "\n";
            fp << "\n";
            fsectionx << "\n";
            fgapcubed << "\n";
            fseparation << n << " " << fCalc.separationIndex() << " "
                        << fCalc.separationX() << " "
                        << fCalc.separationBlendEndX() << " "
                        << fCalc.separationPressure() << "\n";
            fu << t << " "
               << stateL.disp[nearestIdxL].uy - geomL.points[nearestIdxL].y << " "
               << stateR.disp[nearestIdxR].uy - geomR.points[nearestIdxR].y << "\n";
            fdispXY << t << " "
                    << stateL.disp[nearestIdxL].ux - geomL.points[nearestIdxL].x << " "
                    << stateL.disp[nearestIdxL].uy - geomL.points[nearestIdxL].y << " "
                    << stateR.disp[nearestIdxR].ux - geomR.points[nearestIdxR].x << " "
                    << stateR.disp[nearestIdxR].uy - geomR.points[nearestIdxR].y << "\n";
            fpv << n << " " << fCalc.outletPressure() << "\n";
            fuv << n << " " << fCalc.currentUg << "\n";
        }
        if (n % params.nwrite == 0) {
            // q is the committed modal state at t_n, exactly matching the
            // displacement/area/pressure records written above.
            writeModalDiagnostics(n, t, "L", mdataL, stateL, nearestIdxL, surfaceModeRmsL);
            writeModalDiagnostics(n, t, "R", mdataR, stateR, nearestIdxR, surfaceModeRmsR);
        }
        
        int icont_used_this_step = 0;
        bool broke_by_convergence = false;
        bool broke_by_no_contact = false;

        fCalc.resetPreviousContactForce();

        // --- 接触反復計算 ---
        for (int icont = 1; icont <= params.ncont; ++icont) {

            icont_used_this_step = icont;
            total_contact_iterations++;

            // 1. モード力への変換 (L / R)
{            auto t0 = now();
            fCalc.projectLoadsToModes();
            auto t1 = now();
            time_f2mode += elapsed_ms(t0, t1);}
            total_f2mode_calls++;

            if (n % 10 == 0 || t > 0.12) {
                auto [maxFiL, imaxFiL] = maxAbsIndex(fCalc.fiL);
                auto [maxFiR, imaxFiR] = maxAbsIndex(fCalc.fiR);

                fmodedbg << std::scientific << std::setprecision(12)
                        << n << "," << t << "," << icont << ","
                        << "after_f2mode,"
                        << maxFiL << "," << imaxFiL << ","
                        << maxFiR << "," << imaxFiR << ","
                        << fCalc.contactFlag << ","
                        << fCalc.max_force_diff
                        << "\n";
            }

            // Newmark parameters
            const double newmark_beta  = 0.275625;
            const double newmark_gamma = 0.55;

            // 2. 時間積分 (左声帯 L)
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < mdataL.nModes; ++i) {
                double f    = fCalc.fiL[i];                      
                double q    = stateL.q[i];                       
                double qdot = stateL.qdot[i];                    
                double qdd  = stateL.qddot[i];                   
                double omega = omegaL[i];

                double qf, qfdot, qfddot;
                integrator.newmarkStep(f, q, qdot, qdd, params.dt, omega, params.zetaL,
                                       newmark_beta, newmark_gamma, qf, qfdot, qfddot);

                stateL.qf[i]     = qf;
                stateL.qfdot[i]  = qfdot;
                stateL.qfddot[i] = qfddot;
            }
            

            // 3. 時間積分 (右声帯 R)
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < mdataR.nModes; ++i) {
                double f    = fCalc.fiR[i];                      
                double q    = stateR.q[i];                       
                double qdot = stateR.qdot[i];                    
                double qdd  = stateR.qddot[i];                   
                double omega = omegaR[i];

                double qf, qfdot, qfddot;
                integrator.newmarkStep(f, q, qdot, qdd, params.dt, omega, params.zetaR,
                                       newmark_beta, newmark_gamma, qf, qfdot, qfddot);

                stateR.qf[i]     = qf;
                stateR.qfdot[i]  = qfdot;
                stateR.qfddot[i] = qfddot;
            }
            

            // 4. モード変位 → 節点変位 (L / R)
{            auto t0 = now();
            stateL.mode2ufSurface(geomL, mdataL, n+1);
            stateR.mode2ufSurface(geomR, mdataR, n+1);
            auto t1 = now();
            time_mode2uf += elapsed_ms(t0, t1);}
            total_mode2uf_calls += 2;

            // 5. 接触判定とめり込み力計算 (L / R 相対計算)
	{            auto t0 = now();
	            fCalc.applyContactLoads(n, icont);
	            auto t1 = now();
	            time_calcDis += elapsed_ms(t0, t1);}
	            total_calcDis_calls++;
            
            // 収束判定
            if (fCalc.contactFlag
                && fCalc.contact_force_residual < 1.0e-4
                && fCalc.contact_penetration_residual < 1.0e-7) {
                broke_by_convergence = true;
                break; 
            }
            if (!fCalc.contactFlag) {
                broke_by_no_contact = true;
                break;  }
        }

        // applyContactLoads() has just replaced the contact iterate in the
        // surface load.  Solve once more from the committed t_n state so the
        // displacement that is finally committed is produced by that exact
        // final load (also covers contact disappearance and ncont == 0).
{       auto t0 = now();
        fCalc.projectLoadsToModes();
        auto t1 = now();
        time_f2mode += elapsed_ms(t0, t1); }
        ++total_f2mode_calls;

        constexpr double newmark_beta_final = 0.275625;
        constexpr double newmark_gamma_final = 0.55;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < mdataL.nModes; ++i) {
            integrator.newmarkStep(fCalc.fiL[i], stateL.q[i], stateL.qdot[i], stateL.qddot[i],
                                   params.dt, omegaL[i], params.zetaL,
                                   newmark_beta_final, newmark_gamma_final,
                                   stateL.qf[i], stateL.qfdot[i], stateL.qfddot[i]);
        }
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < mdataR.nModes; ++i) {
            integrator.newmarkStep(fCalc.fiR[i], stateR.q[i], stateR.qdot[i], stateR.qddot[i],
                                   params.dt, omegaR[i], params.zetaR,
                                   newmark_beta_final, newmark_gamma_final,
                                   stateR.qf[i], stateR.qfdot[i], stateR.qfddot[i]);
        }
{       auto t0 = now();
        stateL.mode2uf(geomL, mdataL, n + 1);
        stateR.mode2uf(geomR, mdataR, n + 1);
        auto t1 = now();
        time_mode2uf += elapsed_ms(t0, t1); }
        total_mode2uf_calls += 2;

        
        max_icont_used = std::max(max_icont_used, icont_used_this_step);

        if (broke_by_convergence) {
            steps_converged_break++;
        } else if (broke_by_no_contact) {
            steps_no_contact_break++;
        } else {
            steps_reached_ncont++;
        }


        // {   //[DEBUG]
        //     auto [minA, maxA, idxMinA, idxMaxA] = minMaxIndex(fCalc.harea);
        //     auto [maxFiL, imaxFiL] = maxAbsIndex(fCalc.fiL);
        //     auto [maxFiR, imaxFiR] = maxAbsIndex(fCalc.fiR);
        //     double maxQL  = maxAbsVector(stateL.qf);
        //     double maxQdL = maxAbsVector(stateL.qfdot);
        //     double maxQaL = maxAbsVector(stateL.qfddot);
        //     double maxQR  = maxAbsVector(stateR.qf);
        //     double maxQdR = maxAbsVector(stateR.qfdot);
        //     double maxQaR = maxAbsVector(stateR.qfddot);
        //     double maxPredDispL = maxAbsNodeDisp(geomL, stateL);
        //     double maxPredDispR = maxAbsNodeDisp(geomR, stateR);
        //     bool diverged =
        //         !std::isfinite(minA) ||
        //         !std::isfinite(maxA) ||
        //         !std::isfinite(fCalc.currentUg) ||
        //         !std::isfinite(fCalc.currentPg) ||
        //         !std::isfinite(maxAbsPsurf()) ||
        //         !std::isfinite(maxFiL) ||
        //         !std::isfinite(maxFiR) ||
        //         !std::isfinite(maxQL) ||
        //         !std::isfinite(maxQdL) ||
        //         !std::isfinite(maxQaL) ||
        //         !std::isfinite(maxQR) ||
        //         !std::isfinite(maxQdR) ||
        //         !std::isfinite(maxQaR) ||
        //         !std::isfinite(maxPredDispL) ||
        //         !std::isfinite(maxPredDispR);
        //     if (n % 10 == 0 || t > 0.12 || diverged) {
        //         fstepdbg << std::scientific << std::setprecision(12)
        //                 << n << "," << t << ","
        //                 << minA << "," << maxA << "," << idxMinA << ","
        //                 << fCalc.currentUg << "," << fCalc.currentPg << ","
        //                 << (fCalc.Pd.size() > 0 ? fCalc.Pd[0] : 0.0) << ","
        //                 << (fCalc.Pd.size() > 9 ? fCalc.Pd[9] : 0.0) << ","
        //                 << maxAbsPsurf() << ","
        //                 << maxFiL << "," << imaxFiL << ","
        //                 << maxFiR << "," << imaxFiR << ","
        //                 << maxQL << "," << maxQdL << "," << maxQaL << ","
        //                 << maxQR << "," << maxQdR << "," << maxQaR << ","
        //                 << maxPredDispL << "," << maxPredDispR << ","
        //                 << icont_used_this_step << ","
        //                 << fCalc.contactFlag << ","
        //                 << fCalc.max_force_diff << ","
        //                 << diverged
        //                 << "\n";
        //     }
        //     if (diverged) {
        //         std::cerr << "[DIVERGED] step=" << n
        //                 << " t=" << t
        //                 << " minA=" << minA
        //                 << " idxMinA=" << idxMinA
        //                 << " Ug=" << fCalc.currentUg
        //                 << " Pg=" << fCalc.currentPg
        //                 << " maxFiL=" << maxFiL
        //                 << " maxFiR=" << maxFiR
        //                 << " maxQL=" << maxQL
        //                 << " maxQR=" << maxQR
        //                 << std::endl;
        //         break;
        //     }
        // }
    
        // 状態の確定
        stateL.uf2u();
        stateR.uf2u();
        auto t0 = now();

        // 3Dモデル出力
        if (n % 20 == 0) {
            //writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 20);
            num++;
        }

        // The acoustic sample belongs to the same t_n fluid update as the
        // other scalar outputs above.
        soundSignal.push_back(fCalc.outletPressure());
        auto t1 = now();
        time_output += elapsed_ms(t0, t1);

    }
    WavWriter::save(soundSignal, params.dt, (runDir / "test_sound.wav").string());
    
    std::cout << "\n=== Timing Summary ===\n";
    std::cout << "calcArea  : " << time_calcArea  << " ms\n";
    std::cout << "calcForce : " << time_calcForce << " ms\n";
    std::cout << "f2mode    : " << time_f2mode    << " ms\n";
    std::cout << "mode2uf   : " << time_mode2uf   << " ms\n";
    std::cout << "calcDis   : " << time_calcDis   << " ms\n";
    std::cout << "output    : " << time_output    << " ms\n";
    std::cout << "for all   : " << (time_calcArea + time_calcDis + time_calcForce + time_f2mode + time_mode2uf + time_output)/60000 << " min\n"; 

    
    std::cout << "\n=== Contact Iteration Summary ===\n";
    std::cout << "total contact iterations : " << total_contact_iterations << "\n";
    std::cout << "avg icont per step       : "
            << static_cast<double>(total_contact_iterations) / params.nstep << "\n";
    std::cout << "max icont used           : " << max_icont_used << "\n";
    std::cout << "steps no contact break   : " << steps_no_contact_break << "\n";
    std::cout << "steps converged break    : " << steps_converged_break << "\n";
    std::cout << "steps reached ncont      : " << steps_reached_ncont << "\n";
    std::cout << "f2mode calls             : " << total_f2mode_calls << "\n";
    std::cout << "mode2uf calls            : " << total_mode2uf_calls << "\n";
    std::cout << "calcDis calls            : " << total_calcDis_calls << "\n";


    std::cout << "[Simulation] Run complete." << std::endl;
} 

void Simulation::writeVTKCombined(int step, const Geometry& geomL, const State& stateL, 
                                  const Geometry& geomR, const State& stateR, 
                                  const std::string& rdir, int nwrite) {
    // ファイル名 (例: deform_combined0000.vtu)
    std::ostringstream num;
    num << std::setw(4) << std::setfill('0') << step;
    std::string filename = rdir + "/deform_combined" + num.str() + ".vtu";

    std::filesystem::create_directories(rdir);
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    std::cout << "step: " << step * nwrite << std::endl;
    std::cout << "output: " << filename << std::endl;   

    int totalPoints = geomL.nPoints + geomR.nPoints;
    int totalCells = geomL.nCells + geomR.nCells;

    fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    fout << "  <UnstructuredGrid>\n";
    fout << "    <Piece NumberOfPoints=\"" << totalPoints 
         << "\" NumberOfCells=\"" << totalCells << "\">\n";

    // ==========================================
    // 1. Points (頂点座標)
    // ==========================================
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    // 左声帯の座標
    for (int i = 0; i < geomL.nPoints; i++) {
        fout << std::scientific << std::setprecision(15)
             << stateL.disp[i].ux << " " << stateL.disp[i].uy << " " << stateL.disp[i].uz << "\n";
    }
    // 右声帯の座標
    for (int i = 0; i < geomR.nPoints; i++) {
        fout << std::scientific << std::setprecision(15)
             << stateR.disp[i].ux << " " << stateR.disp[i].uy << " " << stateR.disp[i].uz << "\n";
    }
    fout << "        </DataArray>\n";
    fout << "      </Points>\n";

    // ==========================================
    // 2. Cells (セル情報)
    // ==========================================
    fout << "      <Cells>\n";
    
    // --- Connectivity (どの頂点が繋がっているか) ---
    fout << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    // 左声帯
    for (int i = 0; i < geomL.nCells; i++) {
        int nverts = geomL.jtypes[geomL.types[i]];
        for (int j = 0; j < nverts; j++) {
            fout << geomL.connect[i][j] << " ";
        }
        fout << "\n";
    }
    // 右声帯（※左声帯の頂点数 geomL.nPoints 分だけインデックスをズラす）
    for (int i = 0; i < geomR.nCells; i++) {
        int nverts = geomR.jtypes[geomR.types[i]];
        for (int j = 0; j < nverts; j++) {
            fout << geomR.connect[i][j] + geomL.nPoints << " ";
        }
        fout << "\n";
    }
    fout << "        </DataArray>\n";

    // --- Offsets (累計頂点数) ---
    fout << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    long long lastOffsetL = 0;
    // 左声帯
    for (int i = 0; i < geomL.nCells; i++) {
        fout << geomL.offsets[i] << "\n";
        if (i == geomL.nCells - 1) {
            lastOffsetL = geomL.offsets[i]; // 左声帯の最後のオフセット値を記憶
        }
    }
    // 右声帯（※左声帯の最後のオフセット値を足し合わせる）
    for (int i = 0; i < geomR.nCells; i++) {
        fout << lastOffsetL + geomR.offsets[i] << "\n";
    }
    fout << "        </DataArray>\n";

    // --- Types (セルの種類) ---
    fout << "        <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < geomL.nCells; i++) {
        fout << geomL.types[i] << "\n";
    }
    for (int i = 0; i < geomR.nCells; i++) {
        fout << geomR.types[i] << "\n";
    }
    fout << "        </DataArray>\n";
    fout << "      </Cells>\n";

    fout << "    </Piece>\n";
    fout << "  </UnstructuredGrid>\n";
    fout << "</VTKFile>\n";

    fout.close();
}
