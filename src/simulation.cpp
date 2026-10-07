
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
             << "diagnostic_output_interval_steps = "
             << params.diagnosticOutputIntervalSteps << "\n"
             << "perturbation_enabled = " << params.perturbationEnabled << "\n"
             << "perturbation_time_sec = " << params.perturbationTimeSec << "\n"
             << "perturbation_mode_index_1based = " << params.perturbationModeIndex << "\n"
             << "perturbation_amplitude_at_probe_mm = "
             << params.perturbationAmplitudeAtProbeMm << "\n"
             << "perturbation_pattern = " << params.perturbationPattern << "\n"
             << "ap_contact_diagnostic_enabled = " << params.apContactDiagnosticEnabled << "\n"
             << "ap_grid_axis = j indexes ascending clustered physical z levels\n"
             << "ap_contact_analysis_rows_per_end = " << params.apContactAnalysisRowsPerEnd << "\n"
             << "ap_contact_analysis_distance_mm = " << params.apContactAnalysisDistanceMm << "\n"
             << "ap_contact_excluded_rows_per_end = " << params.apContactExcludedRowsPerEnd << "\n"
             << "ap_contact_excluded_distance_mm = " << params.apContactExcludedDistanceMm << "\n"
             << "ap_low_rows_0based = [0," << std::max(0, params.apContactAnalysisRowsPerEnd - 1) << "]\n"
             << "ap_high_rows_0based = [N_AP-" << params.apContactAnalysisRowsPerEnd
             << ",N_AP-1]\n"
             << "contact_detail_start_sec = " << params.contactDetailStartSec << "\n"
             << "contact_detail_end_sec = " << params.contactDetailEndSec << "\n"
             << "diagnostic_load_mode = " << params.diagnosticLoadMode << "\n"
             << "diagnostic_reference_window_start_sec = "
             << params.diagnosticReferenceWindowStartSec << "\n"
             << "diagnostic_reference_window_end_sec = "
             << params.diagnosticReferenceWindowEndSec << "\n"
             << "fixed_node_ids_file = " << params.fixedNodeIdsFile.string() << "\n"
             << "fixed_boundary_diagnostic_scope = "
             << (params.fixedNodeIdsFile.empty()
                    ? "AP end candidates only; no COMSOL fixed IDs supplied"
                    : "explicit zero-based VTU point IDs") << "\n"
             << "left_mode_vtu = " << leftModeFile.string() << "\n"
             << "right_mode_vtu = " << rightModeFile.string() << "\n"
             << "left_frequency = " << leftFrequencyFile.string() << "\n"
             << "right_frequency = " << rightFrequencyFile.string() << "\n"
             << "left_surface_nas = " << leftSurfaceNasFile.string() << "\n"
             << "right_surface_nas = " << rightSurfaceNasFile.string() << "\n"
             << "flow_sections = 50\n"
             << "area_close_m2 = 1e-8\n";
    fCalc.setOutputDirectory(runDir);
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

    // Optional COMSOL/VTU fixed-node audit. IDs are explicitly interpreted as
    // zero-based VTU point IDs and are checked independently on both folds.
    if (!params.fixedNodeIdsFile.empty()) {
        const fs::path fixedPath = resolveInputPath(params.fixedNodeIdsFile);
        requireInputFile(fixedPath, "Fixed-node ID file");
        std::ifstream input(fixedPath);
        std::vector<int> fixedIds;
        std::string line;
        while (std::getline(input, line)) {
            const auto hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream values(line);
            int id = -1;
            while (values >> id) fixedIds.push_back(id);
        }
        if (fixedIds.empty())
            throw std::runtime_error("Fixed-node ID file contains no IDs: " + fixedPath.string());
        std::ofstream audit(runDir / "csv" / "fixed_boundary_mode_audit.csv");
        audit << "side,mode_index_1based,frequency_hz,fixed_node_count,"
              << "max_fixed_mode_displacement,max_all_mode_displacement,fixed_to_all_ratio,id_base\n";
        auto writeAudit = [&](const char* side, const Geometry& geometry,
                              const ModeData& modes) {
            for (int mode = 0; mode < modes.nModes; ++mode) {
                double allMax = 0.0, fixedMax = 0.0;
                int validCount = 0;
                for (const auto& u : modes.modes[mode]) {
                    allMax = std::max(allMax, std::sqrt(u.ux*u.ux + u.uy*u.uy + u.uz*u.uz));
                }
                for (int id : fixedIds) {
                    if (id < 0 || id >= geometry.nPoints) continue;
                    const auto& u = modes.modes[mode][id];
                    fixedMax = std::max(fixedMax, std::sqrt(u.ux*u.ux + u.uy*u.uy + u.uz*u.uz));
                    ++validCount;
                }
                audit << std::scientific << std::setprecision(15)
                      << side << "," << modes.sourceModeIndex(mode) + 1 << ","
                      << modes.frequencies[mode] << "," << validCount << ","
                      << fixedMax << "," << allMax << ","
                      << (allMax > 0.0 ? fixedMax / allMax : 0.0) << ",zero_based_vtu\n";
            }
        };
        writeAudit("L", geomL, mdataL);
        writeAudit("R", geomR, mdataR);
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
    const fs::path csvDir = runDir / "csv";
    std::ofstream fa(datDir / "area.dat");
    std::ofstream fu(datDir / "displace.dat");
    std::ofstream fp(datDir / "pressure.dat");
    std::ofstream fpv(datDir / "pressure_vt.dat");
    std::ofstream fuv(datDir / "airflow_vt.dat");
    std::ofstream fsectionx(datDir / "section_x.dat");
    std::ofstream fgapcubed(datDir / "gap_cubed.dat");
    std::ofstream fseparation(datDir / "separation.dat");
    std::ofstream fdispXY(datDir / "displace_xy.dat");
    std::ofstream fmodal(csvDir / "modal_contribution.csv");
    std::ofstream fmodalDominant(csvDir / "modal_dominant.csv");
    std::ofstream fmodalTop10(csvDir / "modal_top10.csv");
    std::ofstream firregularity(csvDir / "irregularity_timeseries.csv");
    std::ofstream fflowDiagnostics(csvDir / "flow_diagnostics_v2.csv");
    std::ofstream fmodalEnergy(csvDir / "modal_energy.csv");
    std::ofstream fenergySummary(csvDir / "energy_summary.csv");
    std::ofstream fgapField(csvDir / "gap_field.csv");
    fgapField << "snapshot,step,time_s,i_0based,j_0based,x_common_mm,z_common_mm,"
              << "left_y_mm,right_y_mm,signed_gap_y_mm,sampling_method\n";
    std::ofstream fgapFieldSummary(csvDir / "gap_field_summary.csv");
    fgapFieldSummary
        << "snapshot,step,time_s,i_0based,x_common_mm,positive_gap_area_mm2,"
        << "min_signed_gap_mm,negative_gap_sample_fraction,fraction_basis\n";
    std::ofstream fmodalForceSplit(csvDir / "modal_force_split.csv");
    fmodalForceSplit
        << "step,load_time_s,state_time_s,side,mode_index_1based,frequency_hz,q,qdot,qddot,"
        << "fluid_force,contact_force_ap_low,contact_force_ap_high,contact_force_interior,"
        << "total_force,applied_force,structural_restoring_force,structural_damping_force,"
        << "equation_residual,computed_contact_ap_low,computed_contact_ap_high,"
        << "computed_contact_interior,diagnostic_load_mode\n";
    std::ofstream fdiagnosticEvents(csvDir / "diagnostic_events.csv");
    fdiagnosticEvents << "step,time_s,event,detail\n";
    std::ofstream fperturbationEnergy(csvDir / "perturbation_energy.csv");
    fperturbationEnergy
        << "step,load_time_s,state_time_s,side,mode_index_1based,frequency_hz,dq,dv,"
        << "epert_generalized,pfluid_pert_generalized_per_s,"
        << "pcontact_pert_generalized_per_s,pdamp_pert_generalized_per_s,"
        << "wfluid_pert_generalized,wcontact_pert_generalized,wdamp_pert_generalized,"
        << "balance_residual_generalized,equilibrium_residual,event_step_excluded,"
        << "integration_method,units\n";
    std::ofstream fperturbationEnergySummary(csvDir / "perturbation_energy_summary.csv");
    fperturbationEnergySummary
        << "step,state_time_s,side,epert_generalized,wfluid_pert_generalized,"
        << "wcontact_pert_generalized,wdamp_pert_generalized,"
        << "balance_residual_generalized,units\n";
    fflowDiagnostics
        << "step,time_s,input_pressure_pa,ramp_ratio,lung_pressure_pa,"
        << "current_pg_pa,flow_rate_m3_s,dflow_dt_m3_s2,min_area_mm2,"
        << "min_area_index,x_min_area_mm,pressure_at_min_area_pa,"
        << "separation_index,x_separation_mm,pressure_at_separation_pa,"
        << "downstream_pressure_pa,max_abs_surface_pressure_pa,"
        << "max_abs_surface_pressure_index,pressure_recovery_clamp_count,"
        << "has_nonfinite,separation_used_fallback,separation_reason,"
        << "separation_area_ratio,target_separation_area_mm2,actual_separation_area_mm2\n";
    fmodalEnergy
        << "step,time_s,side,mode_index_1based,frequency_hz,q,qdot,"
        << "fluid_modal_force,fluid_power,damping_power,modal_energy,contact_flag\n";
    fenergySummary
        << "step,time_s,fluid_power_left,fluid_power_right,fluid_power_total,"
        << "damping_power_left,damping_power_right,damping_power_total,"
        << "modal_energy_left,modal_energy_right,modal_energy_total,"
        << "cumulative_fluid_work,cumulative_damping_loss,contact_flag\n";
    std::ofstream fcontactCoupling(csvDir / "contact_coupling_debug.csv");
    fcontactCoupling
        << "step,time,contact_iter,max_fiL,max_fiR,"
        << "max_contact_modal_delta_L,max_contact_modal_delta_R,"
        << "max_predicted_xy_change_L_mm,max_predicted_xy_change_R_mm,"
        << "max_pen_m,contact_flag,force_residual,penetration_residual_m\n";
    //[DEBUG]
    std::ofstream fstepdbg(csvDir / "debug_step_summary.csv");
    fstepdbg << "step,time,"
            << "minArea,maxArea,idxMinArea,"
            << "currentUg,currentPg,Pd0,Pd9,"
            << "maxAbsPsurf,"
            << "maxFiL,imaxFiL,maxFiR,imaxFiR,"
            << "maxQL,maxQdL,maxQaL,maxQR,maxQdR,maxQaR,"
            << "maxPredDispL,maxPredDispR,"
            << "icont_used,contactFlag,max_force_diff,"
            << "diverged\n";
    
    std::ofstream fmodedbg(csvDir / "debug_mode_summary.csv");
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
    fseparation << "# step sep_index x_sep[mm] x_blend_end[mm] p_sep[Pa] "
                << "min_area_mm2 target_area_mm2 sep_area_mm2 used_fallback\n";
    fdispXY << "# time[s] uxL[mm] uyL[mm] uxR[mm] uyR[mm]\n";
    fmodal << "step,time_s,side,mode_index,frequency_hz,q,qdot,"
           << "probe_ux_mm,probe_uy_mm,probe_uz_mm,"
           << "surface_rms_ux_mm,surface_rms_uy_mm,surface_rms_uz_mm,"
           << "surface_rms_norm_mm\n";
    fmodalDominant << "step,time_s,side,"
                   << "dominant_probe_uy_mode,dominant_probe_uy_mm,"
                   << "dominant_surface_mode,dominant_surface_rms_mm\n";
    fmodalTop10
        << "step,time_s,side,scope,rank,mode_index,frequency_hz,"
        << "q,modal_norm_mm,magnitude_ratio,signed_projection_ratio,"
        << "cumulative_magnitude_ratio,total_displacement_norm_mm,"
        << "cancellation_factor\n";
    firregularity
        << "step,time_s,x_left_mm,x_right_mm,glottal_area_min_mm2,"
        << "glottal_area_max_mm2,flow_rate_m3_s,outlet_pressure_pa,"
        << "contact_flag\n";


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
    fCalc.setContactMonitor(monitorI, monitorJ);

    int perturbationModeL = -1;
    int perturbationModeR = -1;
    if (params.perturbationEnabled) {
        for (int mode = 0; mode < mdataL.nModes; ++mode) {
            if (mdataL.sourceModeIndex(mode) + 1 == params.perturbationModeIndex) {
                perturbationModeL = mode;
                break;
            }
        }
        for (int mode = 0; mode < mdataR.nModes; ++mode) {
            if (mdataR.sourceModeIndex(mode) + 1 == params.perturbationModeIndex) {
                perturbationModeR = mode;
                break;
            }
        }
        if (perturbationModeL < 0 || perturbationModeR < 0) {
            throw std::runtime_error(
                "perturbationModeIndex is not present in the active mode selection");
        }
    }

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
                   << step << "," << time << "," << side << ","
                   << modes.sourceModeIndex(m) << ","
                   << modes.frequencies[m] << "," << state.q[m] << "," << state.qdot[m] << ","
                   << probeUx << "," << probeUy << "," << probeUz << ","
                   << rmsUx << "," << rmsUy << "," << rmsUz << "," << rmsNorm << "\n";
            if (std::abs(probeUy) > dominantProbeMagnitude) {
                dominantProbeMagnitude = std::abs(probeUy);
                dominantProbeMode = modes.sourceModeIndex(m);
            }
            if (rmsNorm > dominantSurfaceMagnitude) {
                dominantSurfaceMagnitude = rmsNorm;
                dominantSurfaceMode = modes.sourceModeIndex(m);
            }
        }
        fmodalDominant << std::scientific << std::setprecision(12)
                       << step << "," << time << "," << side << ","
                       << dominantProbeMode << "," << dominantProbeMagnitude << ","
                       << dominantSurfaceMode << "," << dominantSurfaceMagnitude << "\n";
    };

    struct ModalContributionEntry {
        int modeIndex = -1;
        double frequencyHz = 0.0;
        double q = 0.0;
        double modalNormMm = 0.0;
        double magnitudeRatio = 0.0;
        double signedProjectionRatio = std::numeric_limits<double>::quiet_NaN();
    };

    auto writeTop10ModalContributions =
        [&](int step, double time, const char* side,
            const ModeData& modes, const State& state, int probeId) {
        constexpr double normEpsMm2 = 1.0e-24;
        constexpr double scalarEpsMm = 1.0e-12;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const int modeCount = modes.nModes;
        const auto& surfacePointIds = state.surfacePointIds;

        std::vector<double> totalUx(surfacePointIds.size(), 0.0);
        std::vector<double> totalUy(surfacePointIds.size(), 0.0);
        std::vector<double> totalUz(surfacePointIds.size(), 0.0);
        for (std::size_t s = 0; s < surfacePointIds.size(); ++s) {
            const int pid = surfacePointIds[s];
            for (int m = 0; m < modeCount; ++m) {
                const double scaleMm = state.q[m] * 1.0e3;
                const auto& phi = modes.modes[m][pid];
                totalUx[s] += scaleMm * phi.ux;
                totalUy[s] += scaleMm * phi.uy;
                totalUz[s] += scaleMm * phi.uz;
            }
        }

        std::vector<double> xyzNormSquared(modeCount, 0.0);
        std::vector<double> xyzProjection(modeCount, 0.0);
        std::vector<double> uyNormSquared(modeCount, 0.0);
        std::vector<double> uyProjection(modeCount, 0.0);
        double totalXyzNormSquared = 0.0;
        double totalUyNormSquared = 0.0;
        for (std::size_t s = 0; s < surfacePointIds.size(); ++s) {
            totalXyzNormSquared += totalUx[s] * totalUx[s]
                                 + totalUy[s] * totalUy[s]
                                 + totalUz[s] * totalUz[s];
            totalUyNormSquared += totalUy[s] * totalUy[s];

            const int pid = surfacePointIds[s];
            for (int m = 0; m < modeCount; ++m) {
                const double scaleMm = state.q[m] * 1.0e3;
                const auto& phi = modes.modes[m][pid];
                const double ux = scaleMm * phi.ux;
                const double uy = scaleMm * phi.uy;
                const double uz = scaleMm * phi.uz;
                xyzNormSquared[m] += ux * ux + uy * uy + uz * uz;
                xyzProjection[m] += ux * totalUx[s]
                                  + uy * totalUy[s]
                                  + uz * totalUz[s];
                uyNormSquared[m] += uy * uy;
                uyProjection[m] += uy * totalUy[s];
            }
        }

        auto makeEntries = [&](const std::vector<double>& normSquared,
                               const std::vector<double>& projection,
                               double totalNormSquared) {
            std::vector<ModalContributionEntry> entries(modeCount);
            double sumModalNorm = 0.0;
            for (int m = 0; m < modeCount; ++m) {
                entries[m].modeIndex = modes.sourceModeIndex(m);
                entries[m].frequencyHz = modes.frequencies[m];
                entries[m].q = state.q[m];
                entries[m].modalNormMm = std::sqrt(normSquared[m]);
                sumModalNorm += entries[m].modalNormMm;
            }
            for (int m = 0; m < modeCount; ++m) {
                if (sumModalNorm > 0.0) {
                    entries[m].magnitudeRatio =
                        entries[m].modalNormMm / sumModalNorm;
                }
                if (totalNormSquared > normEpsMm2) {
                    entries[m].signedProjectionRatio =
                        projection[m] / totalNormSquared;
                }
            }
            return std::make_pair(std::move(entries), sumModalNorm);
        };

        auto writeScope = [&](const char* scope,
                              std::vector<ModalContributionEntry> entries,
                              double totalDisplacementNormMm,
                              double cancellationFactor) {
            std::sort(
                entries.begin(), entries.end(),
                [](const ModalContributionEntry& a,
                   const ModalContributionEntry& b) {
                    if (a.modalNormMm != b.modalNormMm) {
                        return a.modalNormMm > b.modalNormMm;
                    }
                    return a.modeIndex < b.modeIndex;
                });

            const int topCount = std::min<int>(10, entries.size());
            double cumulative = 0.0;
            for (int rank = 0; rank < topCount; ++rank) {
                const auto& entry = entries[rank];
                cumulative += entry.magnitudeRatio;
                fmodalTop10 << std::scientific << std::setprecision(12)
                            << step << "," << time << "," << side << ","
                            << scope << "," << rank + 1 << ","
                            << entry.modeIndex << "," << entry.frequencyHz << ","
                            << entry.q << "," << entry.modalNormMm << ","
                            << entry.magnitudeRatio << ","
                            << entry.signedProjectionRatio << ","
                            << cumulative << "," << totalDisplacementNormMm << ","
                            << cancellationFactor << "\n";
            }
        };

        auto xyzData =
            makeEntries(xyzNormSquared, xyzProjection, totalXyzNormSquared);
        const double totalXyzNormMm = std::sqrt(totalXyzNormSquared);
        const double xyzCancellation =
            totalXyzNormSquared > normEpsMm2
                ? xyzData.second / totalXyzNormMm
                : nan;
        writeScope("surface_xyz", std::move(xyzData.first),
                   totalXyzNormMm, xyzCancellation);

        auto uyData =
            makeEntries(uyNormSquared, uyProjection, totalUyNormSquared);
        const double totalUyNormMm = std::sqrt(totalUyNormSquared);
        const double uyCancellation =
            totalUyNormSquared > normEpsMm2
                ? uyData.second / totalUyNormMm
                : nan;
        writeScope("surface_uy", std::move(uyData.first),
                   totalUyNormMm, uyCancellation);

        std::vector<ModalContributionEntry> probeEntries(modeCount);
        double probeTotalUy = 0.0;
        double sumAbsProbeModeUy = 0.0;
        for (int m = 0; m < modeCount; ++m) {
            const double probeModeUy =
                state.q[m] * 1.0e3 * modes.modes[m][probeId].uy;
            probeEntries[m].modeIndex = modes.sourceModeIndex(m);
            probeEntries[m].frequencyHz = modes.frequencies[m];
            probeEntries[m].q = state.q[m];
            probeEntries[m].modalNormMm = std::abs(probeModeUy);
            probeTotalUy += probeModeUy;
            sumAbsProbeModeUy += std::abs(probeModeUy);
        }
        for (int m = 0; m < modeCount; ++m) {
            if (sumAbsProbeModeUy > 0.0) {
                probeEntries[m].magnitudeRatio =
                    probeEntries[m].modalNormMm / sumAbsProbeModeUy;
            }
            if (std::abs(probeTotalUy) > scalarEpsMm) {
                const double probeModeUy =
                    state.q[m] * 1.0e3 * modes.modes[m][probeId].uy;
                probeEntries[m].signedProjectionRatio =
                    probeModeUy / probeTotalUy;
            }
        }
        const double probeCancellation =
            std::abs(probeTotalUy) > scalarEpsMm
                ? sumAbsProbeModeUy / std::abs(probeTotalUy)
                : nan;
        writeScope("probe_uy", std::move(probeEntries),
                   std::abs(probeTotalUy), probeCancellation);
    };


    stateL.mode2uf(geomL, mdataL, 0); 
    stateL.uf2u(); 
    stateR.mode2uf(geomR, mdataR, 0); 
    stateR.uf2u();

    auto writeGapSnapshot = [&](const char* label, int step, double time) {
        const int ni = std::min(geomL.nxsup, geomR.nxsup);
        const int nj = std::min(geomL.nsurfz, geomR.nsurfz);
        std::vector<std::vector<double>> x(ni, std::vector<double>(nj));
        std::vector<std::vector<double>> z(ni, std::vector<double>(nj));
        std::vector<std::vector<double>> gap(ni, std::vector<double>(nj));
        for (int i = 0; i < ni; ++i) for (int j = 0; j < nj; ++j) {
            const int leftId = geomL.surfp[i][j];
            const int rightId = geomR.surfp[i][j];
            if (leftId < 0 || rightId < 0) continue;
            const auto& left = stateL.disp[leftId];
            const auto& right = stateR.disp[rightId];
            x[i][j] = 0.5 * (left.ux + right.ux);
            z[i][j] = 0.5 * (left.uz + right.uz);
            gap[i][j] = right.uy - left.uy;
            fgapField << std::scientific << std::setprecision(15)
                << label << "," << step << "," << time << "," << i << "," << j << ","
                << x[i][j] << "," << z[i][j] << ","
                << left.uy << "," << right.uy << "," << gap[i][j]
                << ",matched_surface_grid_nodes\n";
        }
        for (int i = 0; i < ni; ++i) {
            double positiveArea = 0.0;
            double minGap = std::numeric_limits<double>::infinity();
            int negative = 0;
            for (int j = 0; j < nj; ++j) {
                minGap = std::min(minGap, gap[i][j]);
                if (gap[i][j] < 0.0) ++negative;
            }
            for (int j = 0; j + 1 < nj; ++j) {
                const double dz = std::abs(z[i][j + 1] - z[i][j]);
                positiveArea += 0.5 * dz
                    * (std::max(0.0, gap[i][j]) + std::max(0.0, gap[i][j + 1]));
            }
            fgapFieldSummary << std::scientific << std::setprecision(15)
                << label << "," << step << "," << time << "," << i << ","
                << x[i][nj / 2] << "," << positiveArea << "," << minGap << ","
                << (nj > 0 ? static_cast<double>(negative) / nj : 0.0)
                << ",matched_grid_node_count\n";
        }
    };
    writeGapSnapshot("initial", 0, 0.0);
    bool referenceGapSnapshotWritten = false;

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

    auto maxAbsDiff = [](const std::vector<double>& a,
                         const std::vector<double>& b) {
        double value = 0.0;
        const std::size_t count = std::min(a.size(), b.size());
        for (std::size_t i = 0; i < count; ++i)
            value = std::max(value, std::abs(a[i] - b[i]));
        return value;
    };

    struct SurfaceXY { double x; double y; };
    auto snapshotSurfaceXY = [](const State& state) {
        std::vector<SurfaceXY> values;
        values.reserve(state.surfacePointIds.size());
        for (int pid : state.surfacePointIds)
            values.push_back({state.predictedDisp[pid].ux, state.predictedDisp[pid].uy});
        return values;
    };
    auto maxSurfaceXYChange = [](const std::vector<SurfaceXY>& a,
                                 const std::vector<SurfaceXY>& b) {
        if (a.size() != b.size() || a.empty())
            return std::numeric_limits<double>::quiet_NaN();
        double value = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const double dx = a[i].x - b[i].x;
            const double dy = a[i].y - b[i].y;
            value = std::max(value, std::sqrt(dx * dx + dy * dy));
        }
        return value;
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

    double cumulativeFluidWork = 0.0;
    double cumulativeDampingLoss = 0.0;
    double previousFluidPowerTotal = 0.0;
    double previousDampingPowerTotal = 0.0;
    bool havePreviousPower = false;
    bool perturbationApplied = false;

    std::vector<double> referenceFluidL(mdataL.nModes, 0.0);
    std::vector<double> referenceFluidR(mdataR.nModes, 0.0);
    std::vector<double> referenceContactL(mdataL.nModes, 0.0);
    std::vector<double> referenceContactR(mdataR.nModes, 0.0);
    std::vector<double> referenceQL(mdataL.nModes, 0.0);
    std::vector<double> referenceQR(mdataR.nModes, 0.0);
    std::vector<double> referenceQdotL(mdataL.nModes, 0.0);
    std::vector<double> referenceQdotR(mdataR.nModes, 0.0);
    std::vector<double> referenceQddotL(mdataL.nModes, 0.0);
    std::vector<double> referenceQddotR(mdataR.nModes, 0.0);
    std::vector<double> cumulativePertFluidL(mdataL.nModes, 0.0);
    std::vector<double> cumulativePertFluidR(mdataR.nModes, 0.0);
    std::vector<double> cumulativePertContactL(mdataL.nModes, 0.0);
    std::vector<double> cumulativePertContactR(mdataR.nModes, 0.0);
    std::vector<double> cumulativePertDampL(mdataL.nModes, 0.0);
    std::vector<double> cumulativePertDampR(mdataR.nModes, 0.0);
    std::vector<double> previousPertFluidPowerL(mdataL.nModes, 0.0);
    std::vector<double> previousPertFluidPowerR(mdataR.nModes, 0.0);
    std::vector<double> previousPertContactPowerL(mdataL.nModes, 0.0);
    std::vector<double> previousPertContactPowerR(mdataR.nModes, 0.0);
    std::vector<double> previousPertDampPowerL(mdataL.nModes, 0.0);
    std::vector<double> previousPertDampPowerR(mdataR.nModes, 0.0);
    std::vector<double> initialPertEnergyL(mdataL.nModes, 0.0);
    std::vector<double> initialPertEnergyR(mdataR.nModes, 0.0);
    long long referenceLoadSamples = 0;
    bool loadSwitchLogged = false;
    bool havePreviousPerturbationPower = false;
    auto applyDiagnosticLoadMode = [&](double time) {
        if (params.diagnosticLoadMode == "live"
            || time < params.perturbationTimeSec) return;
        if (referenceLoadSamples <= 0) {
            throw std::runtime_error(
                "Frozen diagnostic load requested without reference-window samples");
        }
        const double inv = 1.0 / static_cast<double>(referenceLoadSamples);
        for (int mode = 0; mode < mdataL.nModes; ++mode) {
            const double currentContact = fCalc.fiL[mode] - fCalc.fiFluidL[mode];
            const bool freezeFluid = params.diagnosticLoadMode == "freeze_fluid"
                                  || params.diagnosticLoadMode == "freeze_all";
            const bool freezeContact = params.diagnosticLoadMode == "freeze_contact"
                                    || params.diagnosticLoadMode == "freeze_all";
            fCalc.fiL[mode] = (freezeFluid ? referenceFluidL[mode] * inv
                                           : fCalc.fiFluidL[mode])
                            + (freezeContact ? referenceContactL[mode] * inv
                                             : currentContact);
        }
        for (int mode = 0; mode < mdataR.nModes; ++mode) {
            const double currentContact = fCalc.fiR[mode] - fCalc.fiFluidR[mode];
            const bool freezeFluid = params.diagnosticLoadMode == "freeze_fluid"
                                  || params.diagnosticLoadMode == "freeze_all";
            const bool freezeContact = params.diagnosticLoadMode == "freeze_contact"
                                    || params.diagnosticLoadMode == "freeze_all";
            fCalc.fiR[mode] = (freezeFluid ? referenceFluidR[mode] * inv
                                           : fCalc.fiFluidR[mode])
                            + (freezeContact ? referenceContactR[mode] * inv
                                             : currentContact);
        }
    };

    struct StepEnergy {
        double fluidL = 0.0, fluidR = 0.0;
        double dampingL = 0.0, dampingR = 0.0;
        double modalL = 0.0, modalR = 0.0;
    } stepEnergy;

    //writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 1);
    //num++;
    std::cout << "[Simulation] Output step 0 (Initial State)." << std::endl;

    for (int n = 0; n < params.nstep; n++) {
        double t = n * params.dt;

        if (params.perturbationEnabled && !perturbationApplied
            && t >= params.perturbationTimeSec) {
            constexpr double probeModeTolerance = 1.0e-12;
            const double probeUyL =
                mdataL.modes[perturbationModeL][nearestIdxL].uy;
            const double probeUyR =
                mdataR.modes[perturbationModeR][nearestIdxR].uy;
            if (std::abs(probeUyL) < probeModeTolerance
                || std::abs(probeUyR) < probeModeTolerance) {
                throw std::runtime_error(
                    "Cannot apply perturbation: selected mode has near-zero probe uy");
            }
            const double targetLeftMm = -params.perturbationAmplitudeAtProbeMm;
            const double targetRightMm = params.perturbationAmplitudeAtProbeMm;
            const double deltaQL = targetLeftMm / (1.0e3 * probeUyL);
            const double deltaQR = targetRightMm / (1.0e3 * probeUyR);
            stateL.q[perturbationModeL] += deltaQL;
            stateR.q[perturbationModeR] += deltaQR;
            stateL.qf = stateL.q;
            stateR.qf = stateR.q;
            stateL.qfdot = stateL.qdot;
            stateR.qfdot = stateR.qdot;
            stateL.qfddot = stateL.qddot;
            stateR.qfddot = stateR.qddot;
            stateL.mode2uf(geomL, mdataL, n);
            stateR.mode2uf(geomR, mdataR, n);
            stateL.uf2u();
            stateR.uf2u();
            perturbationApplied = true;
            writeGapSnapshot("post_perturbation", n, t);
            fdiagnosticEvents << n << "," << std::scientific << std::setprecision(15)
                              << t << ",perturbation,modal_coordinate_increment\n";

            std::ofstream perturbationManifest(
                runDir / "txt" / "perturbation_manifest.txt", std::ios::trunc);
            perturbationManifest << std::scientific << std::setprecision(15)
                << "requested_time_s = " << params.perturbationTimeSec << "\n"
                << "applied_step = " << n << "\n"
                << "applied_time_s = " << t << "\n"
                << "mode_index_1based = " << params.perturbationModeIndex << "\n"
                << "pattern = " << params.perturbationPattern << "\n"
                << "delta_q_left = " << deltaQL << "\n"
                << "delta_q_right = " << deltaQR << "\n"
                << "realized_probe_increment_left_mm = "
                << 1.0e3 * probeUyL * deltaQL << "\n"
                << "realized_probe_increment_right_mm = "
                << 1.0e3 * probeUyR * deltaQR << "\n";
        }

        if (!referenceGapSnapshotWritten
            && t >= params.diagnosticReferenceWindowEndSec) {
            writeGapSnapshot("reference_window_end", n, t);
            referenceGapSnapshotWritten = true;
        }
        if (!loadSwitchLogged && params.diagnosticLoadMode != "live"
            && t >= params.perturbationTimeSec) {
            fdiagnosticEvents << n << "," << std::scientific << std::setprecision(15)
                              << t << ",load_mode_switch," << params.diagnosticLoadMode << "\n";
            loadSwitchLogged = true;
        }

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

        // Capture the contact-free generalized force before any contact load
        // is assembled.  This projection is diagnostic-only and does not
        // alter the surface load or committed structural state.
        fCalc.projectFluidLoadsToModes();

        stepEnergy = {};
        for (int mode = 0; mode < mdataL.nModes; ++mode) {
            const double q = stateL.q[mode];
            const double qdot = stateL.qdot[mode];
            const double omega = omegaL[mode];
            stepEnergy.fluidL += fCalc.fiFluidL[mode] * qdot;
            stepEnergy.dampingL += 2.0 * params.zetaL * omega * qdot * qdot;
            stepEnergy.modalL += 0.5 * (qdot * qdot + omega * omega * q * q);
        }
        for (int mode = 0; mode < mdataR.nModes; ++mode) {
            const double q = stateR.q[mode];
            const double qdot = stateR.qdot[mode];
            const double omega = omegaR[mode];
            stepEnergy.fluidR += fCalc.fiFluidR[mode] * qdot;
            stepEnergy.dampingR += 2.0 * params.zetaR * omega * qdot * qdot;
            stepEnergy.modalR += 0.5 * (qdot * qdot + omega * omega * q * q);
        }
        const double fluidPowerTotal = stepEnergy.fluidL + stepEnergy.fluidR;
        const double dampingPowerTotal = stepEnergy.dampingL + stepEnergy.dampingR;
        if (havePreviousPower) {
            cumulativeFluidWork += 0.5 * params.dt
                * (previousFluidPowerTotal + fluidPowerTotal);
            cumulativeDampingLoss += 0.5 * params.dt
                * (previousDampingPowerTotal + dampingPowerTotal);
        }
        previousFluidPowerTotal = fluidPowerTotal;
        previousDampingPowerTotal = dampingPowerTotal;
        havePreviousPower = true;

        if (n % params.diagnosticOutputIntervalSteps == 0) {
            int minAreaIndex = -1;
            double minArea = std::numeric_limits<double>::infinity();
            // Match the physical flow constriction search: the first and last
            // stations are boundary planes rather than candidate constrictions.
            for (int index = 1; index + 1 < static_cast<int>(fCalc.harea.size()); ++index) {
                if (fCalc.harea[index] < minArea) {
                    minArea = fCalc.harea[index];
                    minAreaIndex = index;
                }
            }
            int maxPressureIndex = fCalc.psurf.empty() ? -1 : 0;
            double maxAbsPressure = 0.0;
            bool hasNonfinite = fCalc.flowHasNonfinite();
            for (int index = 0; index < static_cast<int>(fCalc.psurf.size()); ++index) {
                if (!std::isfinite(fCalc.psurf[index])) hasNonfinite = true;
                if (std::abs(fCalc.psurf[index]) > maxAbsPressure) {
                    maxAbsPressure = std::abs(fCalc.psurf[index]);
                    maxPressureIndex = index;
                }
            }
            hasNonfinite = hasNonfinite || !std::isfinite(minArea)
                || !std::isfinite(fCalc.sectionX(minAreaIndex))
                || !std::isfinite(fCalc.separationX())
                || !std::isfinite(fCalc.separationPressure())
                || !std::isfinite(fCalc.downstreamPressure());
            fflowDiagnostics << std::scientific << std::setprecision(12)
                << n << "," << t << "," << params.ps << ","
                << fCalc.rampRatio() << "," << fCalc.lungPressure() << ","
                << fCalc.currentPg << "," << fCalc.currentUg << ","
                << (fCalc.currentUg - fCalc.previousUg) / params.dt << ","
                << minArea << "," << minAreaIndex << ","
                << fCalc.sectionX(minAreaIndex) << ","
                << (minAreaIndex >= 0 ? fCalc.psurf[minAreaIndex]
                                      : std::numeric_limits<double>::quiet_NaN()) << ","
                << fCalc.separationIndex() << "," << fCalc.separationX() << ","
                << fCalc.separationPressure() << "," << fCalc.downstreamPressure() << ","
                << maxAbsPressure << "," << maxPressureIndex << ","
                << fCalc.pressureRecoveryClampCount() << ","
                << static_cast<int>(hasNonfinite) << ","
                << static_cast<int>(fCalc.separationUsedFallback()) << ","
                << fCalc.separationReason() << ","
                << params.flowSeparationAreaRatio << ","
                << fCalc.separationTargetArea() << ","
                << fCalc.separationArea() << "\n";
        }

        if (n % 5 == 0) {
            const auto areaExtrema = std::minmax_element(
                fCalc.harea.begin(), fCalc.harea.end());
            const double minGlottalArea =
                areaExtrema.first != fCalc.harea.end()
                    ? *areaExtrema.first
                    : std::numeric_limits<double>::quiet_NaN();
            const double maxGlottalArea =
                areaExtrema.second != fCalc.harea.end()
                    ? *areaExtrema.second
                    : std::numeric_limits<double>::quiet_NaN();
            const double xLeft =
                stateL.disp[nearestIdxL].uy - geomL.points[nearestIdxL].y;
            const double xRight =
                stateR.disp[nearestIdxR].uy - geomR.points[nearestIdxR].y;

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
            firregularity << std::scientific << std::setprecision(12)
                          << n << "," << t << ","
                          << xLeft << "," << xRight << ","
                          << minGlottalArea << "," << maxGlottalArea << ","
                          << fCalc.currentUg << "," << fCalc.outletPressure()
                          << "," << static_cast<int>(fCalc.contactFlag) << "\n";
        }
        if (n % params.nwrite == 0) {
            // q is the committed modal state at t_n, exactly matching the
            // displacement/area/pressure records written above.
            writeModalDiagnostics(n, t, "L", mdataL, stateL, nearestIdxL, surfaceModeRmsL);
            writeModalDiagnostics(n, t, "R", mdataR, stateR, nearestIdxR, surfaceModeRmsR);
            writeTop10ModalContributions(
                n, t, "L", mdataL, stateL, nearestIdxL);
            writeTop10ModalContributions(
                n, t, "R", mdataR, stateR, nearestIdxR);
        }
        
        int icont_used_this_step = 0;
        bool broke_by_convergence = false;
        bool broke_by_no_contact = false;

        fCalc.resetPreviousContactForce();

        std::vector<double> baseFiL, baseFiR;
        std::vector<SurfaceXY> prevPredictedL, prevPredictedR;

        // --- 接触反復計算 ---
        for (int icont = 1; icont <= params.ncont; ++icont) {

            icont_used_this_step = icont;
            total_contact_iterations++;

            // 1. モード力への変換 (L / R)
{            auto t0 = now();
            fCalc.projectLoadsToModes();
            applyDiagnosticLoadMode(t);
            auto t1 = now();
            time_f2mode += elapsed_ms(t0, t1);}
            total_f2mode_calls++;

            if (icont == 1) {
                baseFiL = fCalc.fiL;
                baseFiR = fCalc.fiR;
            }
            const double modalContactDeltaL = baseFiL.empty()
                ? 0.0 : maxAbsDiff(fCalc.fiL, baseFiL);
            const double modalContactDeltaR = baseFiR.empty()
                ? 0.0 : maxAbsDiff(fCalc.fiR, baseFiR);

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

            const auto currentPredictedL = snapshotSurfaceXY(stateL);
            const auto currentPredictedR = snapshotSurfaceXY(stateR);
            const double predictedChangeL = prevPredictedL.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : maxSurfaceXYChange(currentPredictedL, prevPredictedL);
            const double predictedChangeR = prevPredictedR.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : maxSurfaceXYChange(currentPredictedR, prevPredictedR);
            prevPredictedL = currentPredictedL;
            prevPredictedR = currentPredictedR;

            // 5. 接触判定とめり込み力計算 (L / R 相対計算)
	{            auto t0 = now();
	            fCalc.applyContactLoads(n, icont);
	            auto t1 = now();
	            time_calcDis += elapsed_ms(t0, t1);}
	            total_calcDis_calls++;

            fcontactCoupling << std::scientific << std::setprecision(12)
                << n << "," << t << "," << icont << ","
                << maxAbsVector(fCalc.fiL) << "," << maxAbsVector(fCalc.fiR) << ","
                << modalContactDeltaL << "," << modalContactDeltaR << ","
                << predictedChangeL << "," << predictedChangeR << ","
                << fCalc.max_contact_penetration << ","
                << static_cast<int>(fCalc.contactFlag) << ","
                << fCalc.contact_force_residual << ","
                << fCalc.contact_penetration_residual << "\n";
            
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

        const std::vector<double> computedTotalL = fCalc.fiL;
        const std::vector<double> computedTotalR = fCalc.fiR;
        if (t >= params.diagnosticReferenceWindowStartSec
            && t <= params.diagnosticReferenceWindowEndSec) {
            for (int mode = 0; mode < mdataL.nModes; ++mode) {
                referenceFluidL[mode] += fCalc.fiFluidL[mode];
                referenceContactL[mode] += fCalc.fiL[mode] - fCalc.fiFluidL[mode];
                referenceQL[mode] += stateL.q[mode];
                referenceQdotL[mode] += stateL.qdot[mode];
                referenceQddotL[mode] += stateL.qddot[mode];
            }
            for (int mode = 0; mode < mdataR.nModes; ++mode) {
                referenceFluidR[mode] += fCalc.fiFluidR[mode];
                referenceContactR[mode] += fCalc.fiR[mode] - fCalc.fiFluidR[mode];
                referenceQR[mode] += stateR.q[mode];
                referenceQdotR[mode] += stateR.qdot[mode];
                referenceQddotR[mode] += stateR.qddot[mode];
            }
            ++referenceLoadSamples;
        }
        applyDiagnosticLoadMode(t);
        const std::vector<double> appliedTotalL = fCalc.fiL;
        const std::vector<double> appliedTotalR = fCalc.fiR;

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

        if (referenceLoadSamples > 0 && t >= params.perturbationTimeSec) {
            const double invReference = 1.0 / static_cast<double>(referenceLoadSamples);
            const bool eventStep = t <= params.perturbationTimeSec + 0.5 * params.dt;
            struct PerturbationTotals {
                double energy = 0.0, fluidWork = 0.0, contactWork = 0.0,
                       dampingWork = 0.0, balance = 0.0;
            } totalsL, totalsR;
            auto updatePerturbationEnergy = [&]
                (const char* side, const ModeData& modes, const State& state,
                 const std::vector<double>& omega, double zeta,
                 const std::vector<double>& referenceQ,
                 const std::vector<double>& referenceQdot,
                 const std::vector<double>& referenceQddot,
                 const std::vector<double>& referenceFluid,
                 const std::vector<double>& referenceContact,
                 const std::vector<double>& currentFluid,
                 const std::vector<double>& computedTotal,
                 std::vector<double>& cumulativeFluid,
                 std::vector<double>& cumulativeContact,
                 std::vector<double>& cumulativeDamp,
                 std::vector<double>& previousFluidPower,
                 std::vector<double>& previousContactPower,
                 std::vector<double>& previousDampPower,
                 std::vector<double>& initialEnergy,
                 PerturbationTotals& totals) {
                const bool freezeFluid = params.diagnosticLoadMode == "freeze_fluid"
                                      || params.diagnosticLoadMode == "freeze_all";
                const bool freezeContact = params.diagnosticLoadMode == "freeze_contact"
                                        || params.diagnosticLoadMode == "freeze_all";
                for (int mode = 0; mode < modes.nModes; ++mode) {
                    const double qbar = referenceQ[mode] * invReference;
                    const double vbar = referenceQdot[mode] * invReference;
                    const double abar = referenceQddot[mode] * invReference;
                    const double fluidBar = referenceFluid[mode] * invReference;
                    const double contactBar = referenceContact[mode] * invReference;
                    const double computedContact = computedTotal[mode] - currentFluid[mode];
                    const double appliedFluid = freezeFluid ? fluidBar : currentFluid[mode];
                    const double appliedContact = freezeContact ? contactBar : computedContact;
                    const double dq = state.qf[mode] - qbar;
                    const double dv = state.qfdot[mode] - vbar;
                    const double energy = 0.5 * (dv * dv + omega[mode] * omega[mode] * dq * dq);
                    const double fluidPower = (appliedFluid - fluidBar) * dv;
                    const double contactPower = (appliedContact - contactBar) * dv;
                    const double dampingPower = 2.0 * zeta * omega[mode] * dv * dv;
                    if (eventStep) initialEnergy[mode] = energy;
                    if (havePreviousPerturbationPower && !eventStep) {
                        cumulativeFluid[mode] += 0.5 * params.dt
                            * (previousFluidPower[mode] + fluidPower);
                        cumulativeContact[mode] += 0.5 * params.dt
                            * (previousContactPower[mode] + contactPower);
                        cumulativeDamp[mode] += 0.5 * params.dt
                            * (previousDampPower[mode] + dampingPower);
                    }
                    previousFluidPower[mode] = fluidPower;
                    previousContactPower[mode] = contactPower;
                    previousDampPower[mode] = dampingPower;
                    const double balance = energy - initialEnergy[mode]
                        - cumulativeFluid[mode] - cumulativeContact[mode] + cumulativeDamp[mode];
                    const double equilibriumResidual = abar
                        + 2.0 * zeta * omega[mode] * vbar
                        + omega[mode] * omega[mode] * qbar - fluidBar - contactBar;
                    totals.energy += energy;
                    totals.fluidWork += cumulativeFluid[mode];
                    totals.contactWork += cumulativeContact[mode];
                    totals.dampingWork += cumulativeDamp[mode];
                    totals.balance += balance;
                    if (n % params.diagnosticOutputIntervalSteps == 0) {
                        fperturbationEnergy << std::scientific << std::setprecision(15)
                            << n << "," << t << "," << t + params.dt << "," << side << ","
                            << modes.sourceModeIndex(mode) + 1 << "," << modes.frequencies[mode] << ","
                            << dq << "," << dv << "," << energy << "," << fluidPower << ","
                            << contactPower << "," << dampingPower << ","
                            << cumulativeFluid[mode] << "," << cumulativeContact[mode] << ","
                            << cumulativeDamp[mode] << "," << balance << ","
                            << equilibriumResidual << "," << static_cast<int>(eventStep) << ","
                            << "internal_dt_trapezoid_load_tn_state_tn1,generalized_mass_normalized\n";
                    }
                }
            };
            updatePerturbationEnergy("L", mdataL, stateL, omegaL, params.zetaL,
                referenceQL, referenceQdotL, referenceQddotL,
                referenceFluidL, referenceContactL, fCalc.fiFluidL, computedTotalL,
                cumulativePertFluidL, cumulativePertContactL, cumulativePertDampL,
                previousPertFluidPowerL, previousPertContactPowerL, previousPertDampPowerL,
                initialPertEnergyL, totalsL);
            updatePerturbationEnergy("R", mdataR, stateR, omegaR, params.zetaR,
                referenceQR, referenceQdotR, referenceQddotR,
                referenceFluidR, referenceContactR, fCalc.fiFluidR, computedTotalR,
                cumulativePertFluidR, cumulativePertContactR, cumulativePertDampR,
                previousPertFluidPowerR, previousPertContactPowerR, previousPertDampPowerR,
                initialPertEnergyR, totalsR);
            havePreviousPerturbationPower = true;
            if (n % params.diagnosticOutputIntervalSteps == 0) {
                for (const auto& entry : {std::make_pair("L", totalsL),
                                          std::make_pair("R", totalsR)}) {
                    fperturbationEnergySummary << std::scientific << std::setprecision(15)
                        << n << "," << t + params.dt << "," << entry.first << ","
                        << entry.second.energy << "," << entry.second.fluidWork << ","
                        << entry.second.contactWork << "," << entry.second.dampingWork << ","
                        << entry.second.balance << ",generalized_mass_normalized\n";
                }
            }
        }

        if (n % params.diagnosticOutputIntervalSteps == 0) {
            auto writeForceSplit = [&](const char* side, const ModeData& modes,
                                       const State& state, const std::vector<double>& omega,
                                       double zeta, const std::vector<double>& fluid,
                                       const std::vector<double>& low,
                                       const std::vector<double>& high,
                                       const std::vector<double>& interior,
                                       const std::vector<double>& computedLow,
                                       const std::vector<double>& computedHigh,
                                       const std::vector<double>& computedInterior,
                                       const std::vector<double>& computedTotal,
                                       const std::vector<double>& appliedTotal) {
                for (int mode = 0; mode < modes.nModes; ++mode) {
                    const double restoring = omega[mode] * omega[mode] * state.qf[mode];
                    const double damping = 2.0 * zeta * omega[mode] * state.qfdot[mode];
                    const double residual = state.qfddot[mode] + damping + restoring
                                          - appliedTotal[mode];
                    fmodalForceSplit << std::scientific << std::setprecision(15)
                        << n << "," << t << "," << t + params.dt << "," << side << ","
                        << modes.sourceModeIndex(mode) + 1 << "," << modes.frequencies[mode] << ","
                        << state.qf[mode] << "," << state.qfdot[mode] << ","
                        << state.qfddot[mode] << "," << fluid[mode] << ","
                        << low[mode] << "," << high[mode] << "," << interior[mode] << ","
                        << computedTotal[mode] << "," << appliedTotal[mode] << ","
                        << restoring << "," << damping << "," << residual << ","
                        << computedLow[mode] << "," << computedHigh[mode] << ","
                        << computedInterior[mode] << "," << params.diagnosticLoadMode << "\n";
                }
            };
            writeForceSplit("L", mdataL, stateL, omegaL, params.zetaL,
                fCalc.fiFluidL, fCalc.fiContactLowL, fCalc.fiContactHighL,
                fCalc.fiContactInteriorL, fCalc.fiContactComputedLowL,
                fCalc.fiContactComputedHighL, fCalc.fiContactComputedInteriorL,
                computedTotalL, appliedTotalL);
            writeForceSplit("R", mdataR, stateR, omegaR, params.zetaR,
                fCalc.fiFluidR, fCalc.fiContactLowR, fCalc.fiContactHighR,
                fCalc.fiContactInteriorR, fCalc.fiContactComputedLowR,
                fCalc.fiContactComputedHighR, fCalc.fiContactComputedInteriorR,
                computedTotalR, appliedTotalR);
        }

        fCalc.applyContactLoads(n, icont_used_this_step + 1, true);

        if (n % params.diagnosticOutputIntervalSteps == 0) {
            const int contact = static_cast<int>(fCalc.contactFlag);
            for (int mode = 0; mode < mdataL.nModes; ++mode) {
                const double q = stateL.q[mode];
                const double qdot = stateL.qdot[mode];
                const double omega = omegaL[mode];
                const double fluidPower = fCalc.fiFluidL[mode] * qdot;
                const double dampingPower =
                    2.0 * params.zetaL * omega * qdot * qdot;
                const double modalEnergy =
                    0.5 * (qdot * qdot + omega * omega * q * q);
                fmodalEnergy << std::scientific << std::setprecision(12)
                    << n << "," << t << ",L,"
                    << mdataL.sourceModeIndex(mode) + 1 << ","
                    << mdataL.frequencies[mode] << "," << q << "," << qdot << ","
                    << fCalc.fiFluidL[mode] << "," << fluidPower << ","
                    << dampingPower << "," << modalEnergy << "," << contact << "\n";
            }
            for (int mode = 0; mode < mdataR.nModes; ++mode) {
                const double q = stateR.q[mode];
                const double qdot = stateR.qdot[mode];
                const double omega = omegaR[mode];
                const double fluidPower = fCalc.fiFluidR[mode] * qdot;
                const double dampingPower =
                    2.0 * params.zetaR * omega * qdot * qdot;
                const double modalEnergy =
                    0.5 * (qdot * qdot + omega * omega * q * q);
                fmodalEnergy << std::scientific << std::setprecision(12)
                    << n << "," << t << ",R,"
                    << mdataR.sourceModeIndex(mode) + 1 << ","
                    << mdataR.frequencies[mode] << "," << q << "," << qdot << ","
                    << fCalc.fiFluidR[mode] << "," << fluidPower << ","
                    << dampingPower << "," << modalEnergy << "," << contact << "\n";
            }
            fenergySummary << std::scientific << std::setprecision(12)
                << n << "," << t << ","
                << stepEnergy.fluidL << "," << stepEnergy.fluidR << ","
                << stepEnergy.fluidL + stepEnergy.fluidR << ","
                << stepEnergy.dampingL << "," << stepEnergy.dampingR << ","
                << stepEnergy.dampingL + stepEnergy.dampingR << ","
                << stepEnergy.modalL << "," << stepEnergy.modalR << ","
                << stepEnergy.modalL + stepEnergy.modalR << ","
                << cumulativeFluidWork << "," << cumulativeDampingLoss << ","
                << contact << "\n";
        }

        
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
            //writeVTKCombined(num, geomL, stateL, geomR, stateR, vtuResultDir.string(), 20);
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
