
#define _USE_MATH_DEFINES  
#include "Simulation.h"
#include "wavwrite.h"
#include <iostream>
#include <algorithm>
#include <cmath>

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

    geomR.loadFromVTK("../input/M5_test/M5_mode_T3_b16c3.vtu");
    geomR.surfExtractFromNAS("../input/M5_test/M5_surface_T3_d2.nas",69,70);
    geomR.surfArea();

    geomR.surfArea();
    geomR.jtypes[5] = 3;   // 三角形
    geomR.jtypes[9] = 4;   // 四角形
    geomR.jtypes[10] = 4;
    geomR.jtypes[13] = 6;  // 六面体

    mdataR.initialize(params.nmode, geomR);

    mdataR.loadFromVTU("../input/M5_test/M5_mode_T3_b16c3.vtu", geomR);
    mdataR.loadFreqDamping("../input/M5_test/M5_freq_T3_d2_b16c3.txt");

    mdataR.normalizeModes( params.mass, geomR);

    double ymid = geomR.ymid[0]; 

    // ① ジオメトリのYを反転
    for (int i = 0; i < geomR.nPoints; ++i) {
        geomR.points[i].y = 2.0 * ymid - geomR.points[i].y;
    }

    for (int m = 0; m < mdataR.nModes; ++m) {
        for (int i = 0; i < geomR.nPoints; ++i) {
            mdataR.modes[m][i].uy = -mdataR.modes[m][i].uy; 
        }
    }

    stateR.initialize(geomR.nPoints, params.nmode, params.nstep, geomR);

    fCalc.initialize();

    std::ofstream initDbg("../output/debug_step20_30.txt", std::ios::trunc);
    if (initDbg) {
        initDbg << "=== Detailed Parameter Log (Steps 20-30) ===\n";
        initDbg.close();
    }


    std::cout << "[Simulation] Initialization complete." << std::endl;
}

void Simulation::run() {
    std::cout << "[Simulation] Running..." << std::endl;

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
    std::ofstream ff("../output/modeforce.dat");
    std::ofstream fpv("../output/pressure_vt.dat");
    std::ofstream fuv("../output/airflow_vt.dat");

    std::ofstream fdbg("../output/debug_fluid.dat");
    std::ofstream fcOut("../output_force/contact_force_matrix.dat");
    fdbg << "# step time currentUg Pd[0] Pu_last minHarea\n";

    fa << "# x[m]  area[m^2]\n";
    fu << "# x[m]  displaceL displaceR\n";
    fp << "# step  pressure[Pa]\n"; 
    fpv << "# x[m]  pressure[Pa]\n";
    fuv << "# x[m]  airflow[l/s]\n";


    double minDist2 = 1e2;
    int nearestIdxL = -1;
    for (int i = 0; i < geomL.nsurfl; ++i) {
        for (int j = 0; j < geomL.nsurfz; ++j) {
            int idx = geomL.surfp[i][j];
            double dx = geomL.points[idx].x - 10;
            double dz = geomL.points[idx].z - 8.6;
            double dist2 = dx*dx + dz*dz ;


            if (dist2 < minDist2) {
                minDist2 = dist2;
                nearestIdxL = idx;
            }
        }
    }
    int nearestIdxR = nearestIdxL;
    std::cout<<"Monitor Node L idx="<<geomL.points[nearestIdxL].x<<", "<<geomL.points[nearestIdxL].y<<", "<<geomL.points[nearestIdxL].z<<"\n";


    stateL.mode2uf(geomL, mdataL, 0); 
    stateL.uf2u(); 
    stateR.mode2uf(geomR, mdataR, 0); 
    stateR.uf2u();

    

    std::vector<double> soundSignal;
        soundSignal.reserve(params.nstep);

    writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 1);
    num++;
    std::cout << "[Simulation] Output step 0 (Initial State)." << std::endl;

    for (int n = 0; n < params.nstep; n++) {
        double t = n * params.dt;

        // 面積・角度の更新 (左右の相対距離で計算)
        fCalc.calcArea();

        // 圧力の計算と、左右への力の分配
        fCalc.calcForce(t, n);

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

        // --- 接触反復計算 ---
        for (int icont = 1; icont <= params.ncont; ++icont) {

            // 1. モード力への変換 (L / R)
            fCalc.f2mode();

            // Newmark parameters
            const double newmark_beta  = 0.275625;
            const double newmark_gamma = 0.55;

            // 2. 時間積分 (左声帯 L)
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
            stateL.mode2uf(geomL, mdataL, n+1);
            stateR.mode2uf(geomR, mdataR, n+1);

            // 5. 接触判定とめり込み力計算 (L / R 相対計算)
            fCalc.calcDis();
            
            // 収束判定
            if (fCalc.contactFlag && fCalc.max_force_diff < 1.0e-6) { 
                break; 
            }
            if (!fCalc.contactFlag) break;  
        }

        // 変位ログ出力 (左右両方)
        if (n % 5 == 0) {
            fu << n * 1e-5 << " "
               << stateL.predictedDisp[nearestIdxL].uy - geomL.points[nearestIdxL].y << " " // L側変位
               << stateR.predictedDisp[nearestIdxR].uy - geomR.points[nearestIdxR].y << "\n"; // R側変位
        }
    
        // 状態の確定
        stateL.uf2u();
        stateR.uf2u();

        if (n >= 200 && n <= 400) {
            std::ofstream dbgFile("../output/debug_step20_30.txt", std::ios::app);
            if (dbgFile) {
                // 小数点以下12桁まで高精度で出力し、微小なズレを逃さない
                dbgFile << std::scientific << std::setprecision(12);
                dbgFile << "================ Step " << n << " ================\n";
                
                // 1. 流体力学の主要パラメータ
                dbgFile << "[Fluid]\n";
                dbgFile << "  minHarea  : " << fCalc.minHarea[n] << "\n";
                dbgFile << "  currentUg : " << fCalc.currentUg << "\n";
                dbgFile << "  currentPg : " << fCalc.currentPg << "\n";
                // 圧力分布の代表点（メッシュ中央付近）
                int mid_i = geomL.nxsup / 2;
                dbgFile << "  psurf[mid]: " << fCalc.psurf[mid_i] << "\n";

                // 2. モード力と構造力学（1次モード代表）
                dbgFile << "[Structure - Mode 0]\n";
                dbgFile << "  fiL[0]    : " << fCalc.fiL[0] << "\n";
                dbgFile << "  qL[0]     : " << stateL.q[0] << "\n";
                dbgFile << "  qdotL[0]  : " << stateL.qdot[0] << "\n";
                dbgFile << "  qddotL[0] : " << stateL.qddot[0] << "\n";

                // 3. 実際の節点変位（モニター用の nearestIdxL を使用）
                dbgFile << "[Nodal Displacement (nearestIdxL)]\n";
                dbgFile << "  dispL.uy  : " << stateL.disp[nearestIdxL].uy << "\n";
                // 速度や予測変位も確認
                dbgFile << "  velL.uy   : " << stateL.vel[nearestIdxL].uy << "\n";
                
                // 4. 接触とループ制御
                dbgFile << "[System]\n";
                dbgFile << "  contact   : " << (fCalc.contactFlag ? "TRUE" : "FALSE") << "\n";
                
                dbgFile << "---------------------------------------\n";
                dbgFile.close();
            }
        }

        // 3Dモデル出力
        if (n % 20 == 0) {
            //writeVTKCombined(num, geomL, stateL, geomR, stateR, "../result", 20);
            //std::cout << n << "\n";
            //fCalc.outputForceVectors(n);
            num++;
        }

        // 1D流体ログ出力
        if (n % 5 == 0) {
            fpv << std::setw(4) << n << " " << std::setw(8) << fCalc.Pd[9] << "\n";
            fuv << std::setw(4) << n << " " << std::setw(8) << fCalc.currentUg << "\n";

            fdbg << std::setw(6) << n 
                 << " " << std::setw(12) << t 
                 << " " << std::setw(12) << fCalc.currentUg 
                 << " " << std::setw(12) << fCalc.Pd[0]     
                 << " " << std::setw(12) << fCalc.Pu[params.N_sub] 
                 << " " << std::setw(12) << fCalc.minHarea[n] 
                 << "\n";
        }
        soundSignal.push_back(fCalc.Pd[9]);
        
    }

    

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
        fout << std::scientific << std::setprecision(6)
             << stateL.disp[i].ux << " " << stateL.disp[i].uy << " " << stateL.disp[i].uz << "\n";
    }
    // 右声帯の座標
    for (int i = 0; i < geomR.nPoints; i++) {
        fout << std::scientific << std::setprecision(6)
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