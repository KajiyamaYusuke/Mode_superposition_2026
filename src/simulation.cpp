
#define _USE_MATH_DEFINES  
#include "Simulation.h"
#include "wavwrite.h"
#include <iostream>
#include <algorithm>
#include <cmath>

void Simulation::initialize() {
    std::cout << "[Simulation] Initializing..." << std::endl;
    std::string err ="error";

    params.loadFromFile("../input/param.txt", err );

    //geom.loadFromVTK("../input/M5/M5_mode_kawahara_mesh7.vtu");
    geom.loadFromVTK("../input/M5/M5_mode_T2.vtu");
    //geom.loadFromVTK("../input/old/no_mem_mode.vtu");

    //geom.surfExtractFromNAS("/home/kajiyama/code/simulation/input/surface_data_renewal.nas",13,18);
    //geom.surfExtractFromNAS("../input/M5/M5_surface_kawahara_mesh7.nas",64,70);
    geom.surfExtractFromNAS("../input/M5/M5_surface_T2.nas",68,70);
    //geom.surfExtractFromNAS("/home/kajiyama/code/simulation/input/surface_data_old_c.nas",21,30);

    //geom.surfExtract("/home/kajiyama/code/simulation/input/old/surface.txt", 20);
    geom.surfArea();
    geom.print();
 
    geom.jtypes[5] = 3;   // 三角形
    geom.jtypes[9] = 4;   // 四角形
    geom.jtypes[10] = 4;
    geom.jtypes[13] = 6;  // 六面体

    mdata.initialize(params.nmode, geom);

    mdata.loadFromVTU("../input/M5/M5_mode_T2.vtu", geom);
    //mdata.loadFromVTU_old("../input/old/no_mem_mode.vtu", geom);
    mdata.loadFreqDamping("../input/M5/M5_freq_T2.txt");
    //mdata.loadFreqDamping("../input/old/no_mem_frequency.txt");


    mdata.normalizeModes( params.mass, geom);
    

    state.initialize(geom.nPoints, params.nmode, params.nstep, geom);

    // ForceCalculator 初期化
    fCalc.initialize(); 


    std::cout << "[Simulation] Initialization complete." << std::endl;
}

void Simulation::run() {
    std::cout << "[Simulation] Running..." << std::endl;

    // nSteps+1 に対応
    state.qf.resize(mdata.nModes, 0.0);
    state.qfdot.resize(mdata.nModes, 0.0);
    state.qfddot.resize(mdata.nModes, 0.0);

    double P = 1;
    int num = 0;

    std::ofstream fa("../output/area.dat");
    std::ofstream fu("../output/displace.dat");
    std::ofstream fp("../output/pressure.dat");

    fa << "# x[m]  area[m^2]\n";
    fu << "# x[m]  displace\n";
    fp << "# x[m]  pressure[Pa]\n"; 

    std::vector<double> zeta(mdata.nModes, 0);
    double omega1 = 50 * 2 * M_PI;
    double omega2 = 250 * 2 * M_PI;
    double alpha = 2*omega1*omega2*((0.0015*omega2 - 0.01*omega1)/(omega2*omega2 - omega1*omega1));
    double beta = 2*(omega2*0.01 - omega1* 0.0015)/(omega2*omega2 - omega1*omega1);

    double minDist2 = 1e2;
    int nearestIdx = -1;
    for (int i = 0; i < geom.nsurfl; ++i) {
        for (int j = 0; j < geom.nsurfz; ++j) {
            int idx = geom.surfp[i][j];
            double dx = geom.points[idx].x - 10;
            double dz = geom.points[idx].z - 8.6;
            double dist2 = dx*dx + dz*dz ;


            if (dist2 < minDist2) {
                minDist2 = dist2;
                nearestIdx = idx;
            }
        }
    }
    std::cout<<"idx="<<geom.points[nearestIdx].x<<", "<<geom.points[nearestIdx].y<<", "<<geom.points[nearestIdx].z<<"\n";

    for ( int i = 0; i < mdata.nModes; ++i){
        zeta[i] = 1/2*(alpha/(2.0 * M_PI * mdata.frequencies[i]) + beta * 2.0 * M_PI * mdata.frequencies[i]);
    }

    state.mode2uf(geom, mdata, 0); 
    state.uf2u(); // dispを更新
    writeVTK(num, geom, state, "../result", 1); // step 0 を出力
    num++;
    std::cout << "[Simulation] Output step 0 (Initial State). check this if bumpy." << std::endl;

    std::vector<double> soundSignal;
        soundSignal.reserve(params.nstep);

    for (int n = 0
        ; n < params.nstep; n++) {
        double t = n * params.dt;


        // 2. 断面積や角度を更新
        state.calcArea(geom);


        fCalc.calcForce(t, n);

        if (n % 100 == 0) {
            fCalc.outputForceVectors(n);
            for (int i = 0; i<25; ++i){
                //std::cout<<"i = "<<i<<"p = "<<fCalc.psurf[i]<<std::endl;
            }
        }

        fCalc.contactForce();



        if ( n%20 == 0){
            fa <<std::setw(4)<< n;
            fp <<std::setw(4)<< n;
            for (int i = 0; i < geom.nxsup; ++i) {
                
                fa << " " <<std::setw(8)<< state.harea[i] << " ";
                fp << " " <<std::setw(8)<< fCalc.psurf[i] << " ";
            }
            fa << "\n";
            fp << "\n";
            
        }


        for (int icont = 1; icont <= params.ncont; ++icont) {

            // 4. モード力への変換
            fCalc.f2mode();

        // double rampDuration = 0.025; // 0.1秒かけて立ち上げ（状況により0.5など長くする）
        // double rampFactor = 1.0;
        
        // if (t < rampDuration) {
        //     rampFactor = t / rampDuration; 
        //     // 例: t=0なら0倍, t=0.05なら0.5倍, t=0.1以上なら1.0倍
        // }

        // // 計算されたモード力すべてに係数をかける
        // for(int i=0; i<mdata.nModes; ++i) {
        //     fCalc.fi[i] *= rampFactor;
        // }

            // 5. 時間積分（RK4）
/*             for (int i = 0; i < mdata.nModes; i++) {
                double f = fCalc.fi[i];
                double q = state.q[i];
                double qdot = state.qdot[i];
                double omg = 2.0 * M_PI * mdata.frequencies[i];
                double qf, qfdot;

                integrator.rungeStep(f, q, qdot, params.dt, omg, params.zeta, qf, qfdot);

                state.qf[i]    = qf;
                state.qfdot[i] = qfdot;

            } */
            
            // Newmark parameters (average acceleration)
            const double beta  = 0.25;
            const double gamma = 0.5;


            for (int i = 0; i < mdata.nModes; ++i) {
                double f    = fCalc.fi[i];                      // モード力 (tilde f_i)
                double q    = state.q[i];                       // q_n
                double qdot = state.qdot[i];                    // qdot_n
                double qdd  = state.qddot[i];                   // qddot_n (保持しておく)
                double omega = 2.0 * M_PI * mdata.frequencies[i];

                double qf, qfdot, qfddot;
                integrator.newmarkStep(f, q, qdot, qdd, params.dt, omega, params.zeta,
                                    beta, gamma, qf, qfdot, qfddot);

                state.qf[i]    = qf;
                state.qfdot[i] = qfdot;
                state.qfddot[i] = qfddot;
            }




            // 6. モード変位 → 節点変位
            state.mode2uf(geom, mdata, n+1);



            // calculate dissipation force for contact
            fCalc.contactFlag = false;
            fCalc.calcDis();

            if (!fCalc.contactFlag) break;  // contactFlg == false の場合はループを抜ける
        }

        if (n%20 ==0){
            fu << n *1e-5 << " "<<state.predictedDisp[nearestIdx].ufy - geom.points[nearestIdx].y<<" "<<state.predictedDisp[nearestIdx].ufx - geom.points[nearestIdx].x<< "\n";
        }
    
        state.uf2u();

        if( n % 20 == 0){
            writeVTK(num, geom, state, "../result", 200);
            num++;
        }
        soundSignal.push_back(fCalc.currentUg);
    }

    WavWriter::save(soundSignal, params.dt, "../output/output_sound.wav");

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
