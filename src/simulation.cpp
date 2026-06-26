
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
#include <array>
#include <string>
#include <limits>




void writeTraceVTK(const std::vector<std::array<double,3>>& trace,
                   const std::string& filename);
void writePointVTK(double xL, double yL, double zL,
                   double xR, double yR, double zR,
                   int step);
void writeLineVTK(
    const std::vector<std::array<double,3>>& pts0,
    const std::vector<std::array<double,3>>& pts1,
    int step);

Simulation::Simulation()
    : fCalc(geomL, geomR, mdataL, mdataR, stateL, stateR, params)
{}

void Simulation::initialize() {
    std::cout << "[Simulation] Initializing..." << std::endl;
    std::string err ="error";

    params.loadFromFile("../input/param.txt", err );
    params.print();

    geomL.loadFromVTK("../input/M5_test/M5_mode_T3_b8c3.vtu");
    geomL.surfExtractFromNAS("../input/M5_test/M5_surface_T3_d2.nas",69,70);
    geomL.surfArea();

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

    geomR.loadFromVTK("../input/M5_test/M5_mode_T3_b4c3.vtu");
    geomR.surfExtractFromNAS("../input/M5_test/M5_surface_T3_d2.nas",69,70);
    geomR.surfArea();

    geomR.surfArea();

    geomR.jtypes[5] = 3;   // 三角形
    geomR.jtypes[9] = 4;   // 四角形
    geomR.jtypes[10] = 4;
    geomR.jtypes[13] = 6;  // 六面体

    mdataR.initialize(params.nmode, geomR);

    mdataR.loadFromVTU("../input/M5_test/M5_mode_T3_b4c3.vtu", geomR);
    mdataR.loadFreqDamping("../input/M5_test/M5_freq_T3_d2_b4c3.txt");

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

    double P = 1;
    int num = 0;

    std::ofstream fa("../output/area.dat");
    std::ofstream fu("../output/displace.dat");
    std::ofstream fp("../output/pressure.dat");
    std::ofstream fpv("../output/pressure_vt.dat");
    std::ofstream fuv("../output/airflow_vt.dat");

    //std::ofstream fdbg("../output/debug_fluid.dat");
    //std::ofstream contactIterDbg("../output/contact_iteration_debug.csv");

    fa << "# x[m]  area[m^2]\n";
    fu << "# x[m]  displaceL displaceR\n";
    fp << "# step  pressure[Pa]\n"; 
    fpv << "# x[m]  pressure[Pa]\n";
    fuv << "# x[m]  airflow[l/s]\n";


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


    stateL.mode2uf(geomL, mdataL, 0); 
    stateL.uf2u(); 
    stateR.mode2uf(geomR, mdataR, 0); 
    stateR.uf2u();

    fCalc.traceL.clear();
    fCalc.traceR.clear();
    fCalc.traceL.push_back({stateL.disp[nearestIdxL].ux, stateL.disp[nearestIdxL].uy, stateL.disp[nearestIdxL].uz});
    fCalc.traceR.push_back({stateR.disp[nearestIdxR].ux, stateR.disp[nearestIdxR].uy, stateR.disp[nearestIdxR].uz});

    

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

    //writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 1);
    //num++;
    std::cout << "[Simulation] Output step 0 (Initial State)." << std::endl;

    for (int n = 0; n < params.nstep; n++) {
        double t = n * params.dt;

        // 面積・角度の更新 (左右の相対距離で計算)
{        auto t0 = now();
        fCalc.calcArea();
        auto t1 = now();
        time_calcArea += elapsed_ms(t0, t1);}


        // 圧力の計算と、左右への力の分配
{        auto t0 = now();
        fCalc.calcForce(t, n);
        auto t1 = now();
        time_calcForce += elapsed_ms(t0, t1);}

        if (n % 5 == 0) {
            fa << std::setw(4) << n;
            fp << std::setw(4) << n;
            for (int i = 0; i < geomL.nxsup; ++i) {
                fa << " " << std::setw(8) << fCalc.harea[i] << " ";
                fp << " " << std::setw(8) << fCalc.psurf[i] << " ";
            }
            fa << "\n";
            fp << "\n";
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
            fCalc.f2mode();
            auto t1 = now();
            time_f2mode += elapsed_ms(t0, t1);}
            total_f2mode_calls++;

            
            // if (n % 10 == 0) {
            //     std::cout << std::scientific << std::setprecision(12)
            //         << "step=" << n
            //         << " fiL0=" << fCalc.fiL[0]
            //         << " qL0=" << stateL.q[0]
            //         << " qdotL0=" << stateL.qdot[0]
            //         << " powerL0=" << fCalc.fiL[0] * stateL.qdot[0]
            //         << " fiR0=" << fCalc.fiR[0]
            //         << " qR0=" << stateR.q[0]
            //         << " qdotR0=" << stateR.qdot[0]
            //         << " powerR0=" << fCalc.fiR[0] * stateR.qdot[0]
            //         << std::endl;
            // }

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
                double omega = 2.0 * M_PI * mdataL.frequencies[i];

                double qf, qfdot, qfddot;
                integrator.newmarkStep(f, q, qdot, qdd, params.dt, omega, params.zeta,
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
                double omega = 2.0 * M_PI * mdataR.frequencies[i];

                double qf, qfdot, qfddot;
                integrator.newmarkStep(f, q, qdot, qdd, params.dt, omega, params.zeta,
                                       newmark_beta, newmark_gamma, qf, qfdot, qfddot);

                stateR.qf[i]     = qf;
                stateR.qfdot[i]  = qfdot;
                stateR.qfddot[i] = qfddot;
            }
            

            // 4. モード変位 → 節点変位 (L / R)
{            auto t0 = now();
            stateL.mode2uf(geomL, mdataL, n+1);
            stateR.mode2uf(geomR, mdataR, n+1);
            auto t1 = now();
            time_mode2uf += elapsed_ms(t0, t1);}
            total_mode2uf_calls += 2;

            // 5. 接触判定とめり込み力計算 (L / R 相対計算)
	{            auto t0 = now();
	            fCalc.calcDis(n, icont);
	            auto t1 = now();
	            time_calcDis += elapsed_ms(t0, t1);}
	            total_calcDis_calls++;
            
            // 収束判定
            if (fCalc.contactFlag && fCalc.max_force_diff < 1.0e-6) { 
                broke_by_convergence = true;
                break; 
            }
            if (!fCalc.contactFlag) {
                broke_by_no_contact = true;
                break;  }
        }

        
        max_icont_used = std::max(max_icont_used, icont_used_this_step);

        if (broke_by_convergence) {
            steps_converged_break++;
        } else if (broke_by_no_contact) {
            steps_no_contact_break++;
        } else {
            steps_reached_ncont++;
        }


        // 変位ログ出力 (左右両方)
        if (n % 5 == 0) {
            fu << n * params.dt << " "
               << stateL.predictedDisp[nearestIdxL].uy - geomL.points[nearestIdxL].y << " " // L側変位
               << stateR.predictedDisp[nearestIdxR].uy - geomR.points[nearestIdxR].y << "\n"; // R側変位
        }
    
        // 状態の確定
        stateL.uf2u();
        stateR.uf2u();
        // fCalc.traceL.push_back({stateL.disp[nearestIdxL].ux, stateL.disp[nearestIdxL].uy, stateL.disp[nearestIdxL].uz});
        // fCalc.traceR.push_back({stateR.disp[nearestIdxR].ux, stateR.disp[nearestIdxR].uy, stateR.disp[nearestIdxR].uz});

{        auto t0 = now();

        // 3Dモデル出力
        if (n % 20 == 0) {
            //writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 20);
            //writeLineVTK(fCalc.lineStartL, fCalc.lineEndL, num);
            std::cout << n << "\n";
            //fCalc.outputForceVectors(n);
            // writePointVTK(
            //     stateL.disp[nearestIdxL].ux,
            //     stateL.disp[nearestIdxL].uy,
            //     stateL.disp[nearestIdxL].uz,
            //     stateR.disp[nearestIdxR].ux,
            //     stateR.disp[nearestIdxR].uy,
            //     stateR.disp[nearestIdxR].uz,
            //     num
            // );
            num++;
        }

        // 1D流体ログ出力
        if (n % 5 == 0) {
            fpv << std::setw(4) << n << " " << std::setw(8) << fCalc.Pd[9] << "\n";
            fuv << std::setw(4) << n << " " << std::setw(8) << fCalc.currentUg << "\n";
        }
        soundSignal.push_back(fCalc.Pd[9]);
        auto t1 = now();
        time_output += elapsed_ms(t0, t1);}

    }
    //writeTraceVTK(fCalc.traceL, "../output/traceL.vtk");
    //writeTraceVTK(fCalc.traceR, "../output/traceR.vtk");
    WavWriter::save(soundSignal, params.dt, "../output/test_sound.wav");
    
    std::cout << "\n=== Timing Summary ===\n";
    std::cout << "calcArea  : " << time_calcArea  << " ms\n";
    std::cout << "calcForce : " << time_calcForce << " ms\n";
    std::cout << "f2mode    : " << time_f2mode    << " ms\n";
    std::cout << "mode2uf   : " << time_mode2uf   << " ms\n";
    std::cout << "calcDis   : " << time_calcDis   << " ms\n";
    std::cout << "output    : " << time_output    << " ms\n";

    
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

void Simulation::writeVTK(int step, const Geometry& geom, const State& state, const std::string& rdir, int nwrite) {
    // ファイル名
    std::ostringstream num;
    num << std::setw(4) << std::setfill('0') << step;
    std::string filename = rdir + "/deform" + num.str() + ".vtu";

    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    std::cout << "step: " << step * nwrite << std::endl;
    std::cout << "output: " << filename << std::endl;   

    fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    fout << "  <UnstructuredGrid>\n";
    fout << "    <Piece NumberOfPoints=\"" << geom.nPoints 
         << "\" NumberOfCells=\"" << geom.nCells << "\">\n";

    // Points
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nPoints; i++) {
        fout << std::scientific << std::setprecision(6)
             << state.disp[i].ux << " " << state.disp[i].uy << " " << state.disp[i].uz << "\n";
    }
    fout << "        </DataArray>\n";
    fout << "      </Points>\n";

    // Cells
    fout << "      <Cells>\n";
    fout << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nCells; i++) {
        int nverts = geom.jtypes[geom.types[i]];
        for (int j = 0; j < nverts; j++) {
            fout << geom.connect[i][j] << " ";
        }
        fout << "\n";
    }
    fout << "        </DataArray>\n";

    fout << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nCells; i++) {
        fout << geom.offsets[i] << "\n";
    }
    fout << "        </DataArray>\n";

    fout << "        <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nCells; i++) {
        fout << geom.types[i] << "\n";
    }
    fout << "        </DataArray>\n";
    fout << "      </Cells>\n";

    fout << "    </Piece>\n";
    fout << "  </UnstructuredGrid>\n";
    fout << "</VTKFile>\n";

    fout.close();
}

void Simulation::writeVTKCombined(int step, const Geometry& geomL, const State& stateL, 
                                  const Geometry& geomR, const State& stateR, 
                                  const std::string& rdir, int nwrite) {
    // ファイル名 (例: deform_combined0000.vtu)
    std::ostringstream num;
    num << std::setw(4) << std::setfill('0') << step;
    std::string filename = rdir + "/deform_combined" + num.str() + ".vtu";

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

void writeTraceVTK(const std::vector<std::array<double,3>>& trace,
                   const std::string& filename)
{
    std::ofstream fout(filename);

    int n = trace.size();

    fout << "# vtk DataFile Version 3.0\n";
    fout << "Trace\n";
    fout << "ASCII\n";
    fout << "DATASET POLYDATA\n";

    // 頂点
    fout << "POINTS " << n << " float\n";
    for (auto& p : trace) {
        fout << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    if (n >= 2) {
        fout << "LINES 1 " << n + 1 << "\n";
        fout << n << " ";
        for (int i = 0; i < n; ++i) fout << i << " ";
        fout << "\n";
    } else if (n == 1) {
        fout << "VERTICES 1 2\n";
        fout << "1 0\n";
    } else {
        fout << "LINES 0 0\n";
    }

    fout.close();
}

void writePointVTK(double xL, double yL, double zL,
                   double xR, double yR, double zR,
                   int step)
{
    // ファイル名
    std::ostringstream oss;
    oss << "../result/point_" << std::setw(4)
        << std::setfill('0') << step << ".vtk";

    std::ofstream fout(oss.str());

    // ===== ヘッダ =====
    fout << "# vtk DataFile Version 3.0\n";
    fout << "Two Points\n";
    fout << "ASCII\n";
    fout << "DATASET POLYDATA\n";

    // ===== 点データ（左右2点）=====
    fout << "POINTS 2 double\n";
    fout << std::scientific << std::setprecision(15);
    fout << xL << " " << yL << " " << zL << "\n";
    fout << xR << " " << yR << " " << zR << "\n";

    // ===== VERTICES（これが重要）=====
    // 2つの点、それぞれ独立した頂点
    fout << "VERTICES 2 4\n";
    fout << "1 0\n";
    fout << "1 1\n";

    // ===== （任意）ラベル用データ =====
    fout << "POINT_DATA 2\n";
    fout << "SCALARS point_id int 1\n";
    fout << "LOOKUP_TABLE default\n";
    fout << "0\n"; // 左
    fout << "1\n"; // 右

    fout.close();
}

void writeLineVTK(
    const std::vector<std::array<double,3>>& pts0,
    const std::vector<std::array<double,3>>& pts1,
    int step)
{
    // pts0[i] → 始点
    // pts1[i] → 終点（同じiで1本の線）

    std::ostringstream oss;
    oss << "../result/contact_" 
        << std::setw(4) << std::setfill('0') << step 
        << ".vtk";

    std::ofstream fout(oss.str());

    int nLines = pts0.size();
    int nPoints = nLines * 2;

    fout << "# vtk DataFile Version 3.0\n";
    fout << "Contact lines\n";
    fout << "ASCII\n";
    fout << "DATASET POLYDATA\n";

    // ===== POINTS =====
    fout << "POINTS " << nPoints << " double\n";
    fout << std::scientific << std::setprecision(15);

    for (int i = 0; i < nLines; ++i) {
        fout << pts0[i][0] << " " << pts0[i][1] << " " << pts0[i][2] << "\n";
        fout << pts1[i][0] << " " << pts1[i][1] << " " << pts1[i][2] << "\n";
    }

    // ===== LINES =====
    fout << "LINES " << nLines << " " << nLines * 3 << "\n";

    int idx = 0;
    for (int i = 0; i < nLines; ++i) {
        fout << "2 " << idx << " " << idx + 1 << "\n";
        idx += 2;
    }

    fout.close();
}
