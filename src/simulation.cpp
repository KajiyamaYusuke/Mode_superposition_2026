
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
#include <utility>
#include <cstdlib>




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
    // This is a latest-result workspace and is overwritten at each executionh;
    // archival copies are an explicit user action rather than an unbounded
    // accumulation of heavy simulation output.
    const fs::path absoluteParameterFile = fs::absolute(parameterFile).lexically_normal();
    const fs::path projectRoot = absoluteParameterFile.parent_path().parent_path();
    const auto resolveInputPath = [&](const fs::path& configuredPath) {
        if (configuredPath.is_absolute()) return configuredPath.lexically_normal();
        return (absoluteParameterFile.parent_path() / configuredPath).lexically_normal();
    };
    const fs::path leftModeFile = resolveInputPath(params.leftModeFile);
    const fs::path rightModeFile = resolveInputPath(params.rightModeFile);
    const fs::path leftFrequencyFile = resolveInputPath(params.leftFrequencyFile);
    const fs::path rightFrequencyFile = resolveInputPath(params.rightFrequencyFile);
    const fs::path leftSurfaceNasFile = resolveInputPath(params.leftSurfaceNasFile);
    const fs::path rightSurfaceNasFile = resolveInputPath(params.rightSurfaceNasFile);
    const auto requireInputFile = [](const fs::path& path, const char* label) {
        if (!fs::is_regular_file(path)) {
            throw std::runtime_error(std::string(label) + " does not exist: " + path.string());
        }
    };
    requireInputFile(leftModeFile, "Left mode file");
    requireInputFile(rightModeFile, "Right mode file");
    requireInputFile(leftFrequencyFile, "Left frequency file");
    requireInputFile(rightFrequencyFile, "Right frequency file");
    requireInputFile(leftSurfaceNasFile, "Left surface NAS file");
    requireInputFile(rightSurfaceNasFile, "Right surface NAS file");
    const double targetInitialGapMm = params.initialGapMm;
    const char* configuredRunDir = std::getenv("SIMULATION_RUN_DIR");
    runDir = configuredRunDir && *configuredRunDir
        ? fs::absolute(configuredRunDir).lexically_normal()
        : projectRoot / "output";
    vtuResultDir = configuredRunDir && *configuredRunDir
        ? runDir / "result" : projectRoot / "result";
    fs::create_directories(runDir);
    for (const char* subdir : {"csv", "dat", "json", "png", "txt", "wav"}) {
        fs::create_directories(runDir / subdir);
    }
    fs::copy_file(absoluteParameterFile, runDir / "txt" / "params_used.txt",
                  fs::copy_options::overwrite_existing);
    std::ofstream manifest(runDir / "txt" / "manifest.txt");
    manifest << "parameter_file = " << absoluteParameterFile << "\n"
             << "nstep = " << params.nstep << "\n"
             << "dt_s = " << params.dt << "\n"
             << "nmode = " << params.nmode << "\n"
             << "mode_selection = ";
    if (params.modeSelection.empty()) {
        manifest << "all (1..nmode)";
    } else {
        for (std::size_t i = 0; i < params.modeSelection.size(); ++i) {
            if (i) manifest << ",";
            manifest << params.modeSelection[i];
        }
    }
    manifest << "\n"
             << "zetaL = " << params.zetaL << "\n"
             << "zetaR = " << params.zetaR << "\n"
             << "kc1 = " << params.kc1 << "\n"
             << "kc2 = " << params.kc2 << "\n"
             << "kc3 = " << params.kc3 << "\n"
             << "force_frequency_hz = " << params.forcef << "\n"
             << "force_amplitude = " << params.famp << "\n"
             << "force_direction = " << params.forceDirection << "\n"
             << "ncont = " << params.ncont << "\n"
             << "iforce = " << params.iforce << "\n"
             << "contact_reference_frequency_hz = " << params.contactReferenceFrequencyHz << "\n"
             << "flow_blend_length_mm = " << params.flowBlendLengthMm << "\n"
             << "flow_separation_method = downstream_area_ratio\n"
             << "flow_separation_area_ratio = " << params.flowSeparationAreaRatio << "\n"
             << "target_initial_minimum_gap_mm = " << targetInitialGapMm << "\n"
             << "initial_gap_note = "
             << (targetInitialGapMm < 0.0
                    ? "diagnostic geometric precompression; not a physical adduction force"
                    : "geometric minimum-gap target") << "\n"
             << "pressure_ramp_time_sec = "
             << (params.pressureRampTimeSecExplicit ? params.pressureRampTimeSec : 0.05) << "\n"
             << "pressure_ramp_option_explicit = " << params.pressureRampTimeSecExplicit << "\n"
             << "ap_contact_diagnostic_enabled = " << params.apContactDiagnosticEnabled << "\n"
             << "ap_grid_axis = j indexes ascending clustered physical z levels\n"
             << "ap_contact_excluded_rows_per_end = " << params.apContactExcludedRowsPerEnd << "\n"
             << "ap_contact_excluded_distance_mm = " << params.apContactExcludedDistanceMm << "\n"
             << "left_mode_vtu = " << leftModeFile.string() << "\n"
             << "right_mode_vtu = " << rightModeFile.string() << "\n"
             << "left_frequency = " << leftFrequencyFile.string() << "\n"
             << "right_frequency = " << rightFrequencyFile.string() << "\n"
             << "left_surface_nas = " << leftSurfaceNasFile.string() << "\n"
             << "right_surface_nas = " << rightSurfaceNasFile.string() << "\n"
             << "flow_sections = 50\n"
             << "area_close_m2 = 1e-8\n";
    std::cout << "[Simulation] Latest-result directory: " << runDir << "\n";

    geomL.loadFromVTK(leftModeFile.string());
    geomL.surfExtractFromNAS(leftSurfaceNasFile.string(),64,69);
    geomL.surfArea();
    geomL.print();
    geomL.jtypes[5] = 3;   // 三角形
    geomL.jtypes[9] = 4;   // 四角形
    geomL.jtypes[10] = 4;
    geomL.jtypes[13] = 6;  // 六面体

    mdataL.initialize(params.nmode, geomL, params.modeSelection);

    mdataL.loadFromVTU(leftModeFile.string(), geomL);
    mdataL.loadFreqDamping(leftFrequencyFile.string());

    mdataL.normalizeModes( params.mass, geomL);

    geomR.loadFromVTK(rightModeFile.string());
    geomR.surfExtractFromNAS(rightSurfaceNasFile.string(),64,69);
    geomR.surfArea();

    geomR.jtypes[5] = 3;   // 三角形
    geomR.jtypes[9] = 4;   // 四角形
    geomR.jtypes[10] = 4;
    geomR.jtypes[13] = 6;  // 六面体

    mdataR.initialize(params.nmode, geomR, params.modeSelection);

    mdataR.loadFromVTU(rightModeFile.string(), geomR);
    mdataR.loadFreqDamping(rightFrequencyFile.string());

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

    // Set the minimum initial medial gap by translating both complete folds
    // symmetrically. Modal vectors are unchanged by this rigid translation.
    double meanGap = 0.0;
    int gapCount = 0;
    const int commonI = std::min(geomL.nxsup, geomR.nxsup);
    const int commonJ = std::min(geomL.nsurfz, geomR.nsurfz);
    for (int i = 0; i < commonI; ++i) {
        for (int j = 0; j < commonJ; ++j) {
            const int leftId = geomL.surfp[i][j];
            const int rightId = geomR.surfp[i][j];
            if (leftId < 0 || rightId < 0) continue;
            meanGap += geomR.points[rightId].y - geomL.points[leftId].y;
            ++gapCount;
        }
    }
    if (gapCount == 0) {
        throw std::runtime_error("Cannot set initial gap: no paired surface points");
    }
    const double gapSign = meanGap < 0.0 ? -1.0 : 1.0;
    double minimumGapMm = std::numeric_limits<double>::infinity();
    for (int i = 0; i < commonI; ++i) {
        for (int j = 0; j < commonJ; ++j) {
            const int leftId = geomL.surfp[i][j];
            const int rightId = geomR.surfp[i][j];
            if (leftId < 0 || rightId < 0) continue;
            minimumGapMm = std::min(
                minimumGapMm,
                gapSign * (geomR.points[rightId].y - geomL.points[leftId].y));
        }
    }
    const double gapCorrectionMm = targetInitialGapMm - minimumGapMm;
    for (auto& point : geomL.points) {
        point.y -= 0.5 * gapSign * gapCorrectionMm;
    }
    for (auto& point : geomR.points) {
        point.y += 0.5 * gapSign * gapCorrectionMm;
    }
    std::cout << "[Simulation] Initial minimum gap: " << minimumGapMm
              << " mm -> " << targetInitialGapMm << " mm\n";

    stateL.initialize(geomL.nPoints, mdataL.nModes, params.nstep, geomL);
    stateR.initialize(geomR.nPoints, mdataR.nModes, params.nstep, geomR);

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

    const fs::path datDir = runDir / "dat";
    std::ofstream fa(datDir / "area.dat");
    std::ofstream fu(datDir / "displace.dat");
    std::ofstream fp(datDir / "pressure.dat");
    std::ofstream fpv(datDir / "pressure_vt.dat");
    std::ofstream fuv(datDir / "airflow_vt.dat");
    std::ofstream fsectionx(datDir / "section_x.dat");
    std::ofstream fgapcubed(datDir / "gap_cubed.dat");
    std::ofstream fseparation(datDir / "separation.dat");
    std::ofstream fdispXY(datDir / "displace_xy.dat");
    fa << "# step  area[mm^2]\n";
    fu << "# time_s  uyL_mm  uyR_mm\n";
    fp << "# step  pressure[Pa]\n"; 
    fpv << "# step  outlet_pressure[Pa]\n";
    fuv << "# step  airflow[m^3/s]\n";
    fsectionx << "# step  section_x[mm]\n";
    fgapcubed << "# step  integral_g_positive_cubed[mm^4]\n";
    fseparation << "# step sep_index x_sep[mm] x_blend_end[mm] p_sep[Pa] "
                << "min_area_mm2 target_area_mm2 sep_area_mm2 used_fallback\n";
    fdispXY << "# time[s] uxL[mm] uyL[mm] uxR[mm] uyR[mm]\n";
    const double monitorTargetX = 9.2;
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

    stateL.mode2uf(geomL, mdataL, 0); 
    stateL.uf2u(); 
    stateR.mode2uf(geomR, mdataR, 0); 
    stateR.uf2u();

    std::vector<double> soundSignal;
        soundSignal.reserve(params.nstep);

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
                        << fCalc.separationPressure() << " "
                        << fCalc.separationMinArea() << " "
                        << fCalc.separationTargetArea() << " "
                        << fCalc.separationArea() << " "
                        << static_cast<int>(fCalc.separationUsedFallback()) << "\n";
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

        fCalc.applyContactLoads(n, icont_used_this_step + 1, true);


        
        max_icont_used = std::max(max_icont_used, icont_used_this_step);

        if (broke_by_convergence) {
            steps_converged_break++;
        } else if (broke_by_no_contact) {
            steps_no_contact_break++;
        } else {
            steps_reached_ncont++;
        }
    
        // 状態の確定
        stateL.uf2u();
        stateR.uf2u();
        auto t0 = now();

        // 3Dモデル出力
        if (n % 20 == 0 && params.nstep-n <= 5000) {
            writeVTKCombined(num, geomL, stateL, geomR, stateR, vtuResultDir.string(), 20);
            num++;
        }

        // The acoustic sample belongs to the same t_n fluid update as the
        // other scalar outputs above.
        soundSignal.push_back(fCalc.outletPressure());
        auto t1 = now();
        time_output += elapsed_ms(t0, t1);

    }
    WavWriter::save(soundSignal, params.dt, (runDir / "wav" / "test_sound.wav").string());
    
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
