#include "ForceCalculator.h"
#include "Displacement.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>

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
      stateL(stL_), stateR(stR_), sp(sp_) {}



void ForceCalculator::initialize() {
    // --- 左右それぞれのサイズを取得 ---
    int nsurfl_L = geomL.nsurfl;
    int nsurfz_L = geomL.nsurfz;
    int nModes_L = modeDataL.nModes;

    int nsurfl_R = geomR.nsurfl;
    int nsurfz_R = geomR.nsurfz;
    int nModes_R = modeDataR.nModes;

    // --- 左声帯 (L) の変数初期化 ---
    fxL.assign(nsurfl_L, std::vector<double>(nsurfz_L, 0.0));
    fyL.assign(nsurfl_L, std::vector<double>(nsurfz_L, 0.0));
    fzL.assign(nsurfl_L, std::vector<double>(nsurfz_L, 0.0));
    fiL.assign(nModes_L, 0.0);

    // --- 右声帯 (R) の変数初期化 ---
    fxR.assign(nsurfl_R, std::vector<double>(nsurfz_R, 0.0));
    fyR.assign(nsurfl_R, std::vector<double>(nsurfz_R, 0.0));
    fzR.assign(nsurfl_R, std::vector<double>(nsurfz_R, 0.0));
    fiR.assign(nModes_R, 0.0);

    // --- 内部バッファ (流路断面数は左右で共通と仮定) ---
    nxsup = geomL.nxsup; 
    fdisL.assign(nxsup, std::vector<double>(geomL.nsurfz - 1, 0.0));
    fdisR.assign(nxsup, std::vector<double>(geomR.nsurfz - 1, 0.0));

    // --- 共有変数の初期化 ---
    psurf.assign(nxsup, 0.0);
    Ug.assign(stateL.nSteps, 0.0);       // ステップ数はL/R共通
    minHarea.assign(stateL.nSteps, 0.0);

    int Nsecg = sp.N_sub;
    int Nsecp = sp.N_vt;

    Uu.assign(Nsecg + 1, 0.0); // +2 for boundaries
    Pu.assign(Nsecg + 2, 0.0);
    Ud.assign(Nsecp, 0.0);
    Pd.assign(Nsecp , 0.0);

    harea.clear();
    degreeL.clear();
    degreeR.clear();
    
    currentUg = 0.0;    

    rho = sp.rho ;
    mu  = sp.mu ;
    lg  = geomL.zmax ; // 声門長は左を基準（左右同じはず）

    c_sound = sp.c_sound;

    double A_inlet = M_PI * sp.r_inlet * sp.r_inlet;
    double A_sub   = M_PI * sp.r_sub * sp.r_sub;
    double A_vt    = M_PI * sp.r_vt * sp.r_vt;

    const double A_vt_thresh = 1e-12;
    hasVocalTract = (sp.L_vt > 1e-6) && (A_vt > A_vt_thresh);
    
    // --- インピーダンス計算 ---
    Lui = rho * sp.L_inlet / (2.0 * A_inlet); 
    Cui = sp.L_inlet * A_inlet / (rho * c_sound * c_sound);

    double dx_sub = sp.L_sub / std::max(1, Nsecg);
    Lu = rho * dx_sub / (2.0 * A_sub); // 元コードの /2.0 を再現
    Cu = dx_sub * A_sub / (rho * c_sound * c_sound);

    alpha1 = -2.5e-5*sp.ps+0.185;
    alpha2 = 1.6e-3*sp.ps+0.6;
    beta = 1.125e-4 * sp.ps + 0.1375;

    // safe R2
    if (A_vt > A_vt_thresh) {
        R2 = alpha1 / (A_sub * A_sub) * std::sqrt(rho * mu * c_sound);
    } else {
        R2 = 0.0;
    }

    if (hasVocalTract) {
        double dx_vt = sp.L_vt / std::max(1, Nsecp); 
        La = rho * dx_vt / A_vt;
        Ca = dx_vt * A_vt / (rho * c_sound * c_sound);
        
        Lr = rho * 1.1 * std::sqrt(A_vt / M_PI) / A_vt;
        Rr = alpha2 * rho * c_sound / (9.0 * M_PI * M_PI * A_vt);
    } else {
        La = 1e-1; // ダミー
        Ca = 1.0e30; 
        Lr = 0.0;
        Rr = 0.0;
    }

    contactForceL_ij.resize(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0));
    contactForceR_ij.resize(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0));

    harea.assign(nxsup, 0.0);
    degreeL.assign(2, std::vector<std::vector<double>>(geomL.nxsup, std::vector<double>(geomL.nsurfz, 0.0)));
    degreeR.assign(2, std::vector<std::vector<double>>(geomR.nxsup, std::vector<double>(geomR.nsurfz, 0.0)));

    debugForceFile.open("../output/debug_force.txt", std::ios::trunc);
    if (debugForceFile) {
        debugForceFile << "=== Debug Force Log Initialized ===\n";
    }
    
    std::cout << "[ForceCalculator] initialized (Asymmetric): "
              << "nModesL=" << nModes_L << ", nModesR=" << nModes_R
              << ", nxsup=" << nxsup 
              << ", L_sub=" << sp.L_sub
              << ", L_vt=" << sp.L_vt
              << std::endl;
}

void ForceCalculator::calcForce(double t, int n) {
    // 左右で共通の分割数を使用
    int nsurfl = geomL.nsurfl;
    int nsurfz = geomL.nsurfz;
    int nxsup  = geomL.nxsup;

    for (int i = 0; i < nsurfl; ++i) {
        std::fill(fxL[i].begin(), fxL[i].end(), 0.0);
        std::fill(fyL[i].begin(), fyL[i].end(), 0.0);
        std::fill(fzL[i].begin(), fzL[i].end(), 0.0);

        std::fill(fxR[i].begin(), fxR[i].end(), 0.0);
        std::fill(fyR[i].begin(), fyR[i].end(), 0.0);
        std::fill(fzR[i].begin(), fzR[i].end(), 0.0);
    }

    for (int i = 0; i < nxsup; i++) {
        for (int j = 0; j < nsurfz - 1; j++) {
            fdisL[i][j] = 0.0;
            fdisR[i][j] = 0.0;
        }
    }

    if (sp.iforce == 1) {
        // ==== sin波加振 (テスト用) ====
        minHarea[n] = *std::min_element(harea.begin(), harea.end());
        
        for (int i = 1; i < nxsup-1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {
                fxL[i][j] = sp.famp * std::sin(2.0 * M_PI * sp.forcef * t);
                fxR[i][j] = sp.famp * std::sin(2.0 * M_PI * sp.forcef * t);
            }
        }
    } else if (sp.iforce == 0) {
        // ==== 1D flow model ====
        double minA = findMinHarea();
        minHarea[n] = minA;


        previousUg = (n > 0 && n-1 < (int)Ug.size()) ? Ug[n-1] : 0.0;

        // 安全な面積を使って流体計算を進める
        calcFlowStep(t, sp.dt, minA * 1e-6);
        
        // 剥離点の特定
        int nsep = findNsep(minA);

        Ug[n] = currentUg;

        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = currentPg;

        if(debugForceFile){
            debugForceFile << "Step: " << std::setw(4) << n 
                           << " | minA: " << std::scientific << std::setprecision(12) << minA 
                           << " | subPg: " << std::scientific << std::setprecision(12) << currentPg  << "\n";
        }

        if (minA > 1e-6) {
            for (int i = 1; i < nxsup; i++) {
                // dxの計算は左声帯の座標を使用
                double dx = std::abs(geomL.points[geomL.surfp[i][int(nsurfz/2)]].x - geomL.points[geomL.surfp[i-1][int(nsurfz/2)]].x);
                double h  = (harea[i] + harea[i-1]) / (2.0 * lg);
                double h_prev = std::max(harea[i-1], 1e-3);
                double h_curr = std::max(harea[i], 1e-3);

                double Ugm = currentUg * 1e6;

                double bernoulli_term = 0.5 * rho * Ugm * Ugm * (1.0/(h_prev*h_prev) - 1.0/(h_curr*h_curr));
                double viscous_term   = 12.0 * mu * dx * Ugm / (lg * pow(h, 3.0)) * 1e3;

                psurf[i] = psurf[i-1] + bernoulli_term - viscous_term;

                if (psurf[i] > psurf[i-1]) { break; } // 圧力回復で剥離
            }
        } else {
            for (int i = 1; i < nsep-1; ++i) psurf[i] = currentPg;
            for (int i = nsep-1; i < nxsup; ++i) psurf[i] = Pd[0];
        }

        for (int i = 1; i < nsep-1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {
                
                int pid_ip1_L = geomL.surfp[i+1][j];
                int pid_im1_L = geomL.surfp[i-1][j];
                int pid_jp1_L = geomL.surfp[i][j+1];
                int pid_jm1_L = geomL.surfp[i][j-1];

                if (pid_ip1_L >= 0 && pid_im1_L >= 0 && pid_jp1_L >= 0 && pid_jm1_L >= 0) {
                    double dxL = 0.5 * (stateL.disp[pid_ip1_L].ux - stateL.disp[pid_im1_L].ux);
                    double dyL = 0.5 * (stateL.disp[pid_ip1_L].uy - stateL.disp[pid_im1_L].uy);
                    double dsL = std::sqrt(dxL*dxL + dyL*dyL);
                    double dzL = 0.5 * (stateL.disp[pid_jp1_L].uz - stateL.disp[pid_jm1_L].uz);

                    fxL[i][j] = psurf[i] * dsL * dzL * 1.0e-6 * std::cos(degreeL[1][i][j]) * std::sin(degreeL[0][i][j]);

                    fyL[i][j] = -psurf[i] * dsL * dzL * 1.0e-6 * std::cos(degreeL[1][i][j]) * std::cos(degreeL[0][i][j]);
                    fzL[i][j] = psurf[i] * dsL * dzL * 1.0e-6 * std::sin(degreeL[1][i][j]);

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
                    double dxR = 0.5 * (stateR.disp[pid_ip1_R].ux - stateR.disp[pid_im1_R].ux);
                    double dyR = 0.5 * (stateR.disp[pid_ip1_R].uy - stateR.disp[pid_im1_R].uy);
                    double dsR = std::sqrt(dxR*dxR + dyR*dyR);
                    double dzR = 0.5 * (stateR.disp[pid_jp1_R].uz - stateR.disp[pid_jm1_R].uz);

                    fxR[i][j] = -psurf[i] * dsR * dzR * 1.0e-6 * std::cos(degreeR[1][i][j]) * std::sin(degreeR[0][i][j]);

                    fyR[i][j] = psurf[i] * dsR * dzR * 1.0e-6 * std::cos(degreeR[1][i][j]) * std::cos(degreeR[0][i][j]);
                    fzR[i][j] = -psurf[i] * dsR * dzR * 1.0e-6 * std::sin(degreeR[1][i][j]);

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

    }
}

void ForceCalculator::f2mode() {

    for (int imode = 0; imode < modeDataL.nModes; imode++) {
        double force = 0.0;
        const auto& mode = modeDataL.modes[imode];
        for (int i = 0; i < geomL.nsurfl; i++) {
            for (int j = 0; j < geomL.nsurfz; j++) {
                int pid = geomL.surfp[i][j];
                if (pid < 0) continue; // 節点が存在しない場合スキップ

                const auto& modeAtPoint = mode[pid];
                force += fxL[i][j] * modeAtPoint.ux
                       + fyL[i][j] * modeAtPoint.uy
                       + fzL[i][j] * modeAtPoint.uz;
            }
        }
        fiL[imode] = force;
    }

    for (int imode = 0; imode < modeDataR.nModes; imode++) {
        double force = 0.0;
        const auto& mode = modeDataR.modes[imode];
        for (int i = 0; i < geomR.nsurfl; i++) {
            for (int j = 0; j < geomR.nsurfz; j++) {
                int pid = geomR.surfp[i][j];
                if (pid < 0) continue; // 節点が存在しない場合スキップ

                const auto& modeAtPoint = mode[pid];
                force += fxR[i][j] * modeAtPoint.ux
                       + fyR[i][j] * modeAtPoint.uy
                       + fzR[i][j] * modeAtPoint.uz;
            }
        }
        fiR[imode] = force;
    }
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

void ForceCalculator::calcArea() {
    // 左右どちらかのメッシュサイズを基準とする（今回は左を基準）
    if (geomL.nxsup < 2 || geomL.nsurfz < 2) return;

    std::fill(harea.begin(), harea.end(), 0.0);
    
    for (int i = 0; i < geomL.nxsup; ++i) {
        double hi = 0.0;
        for (int j = 0; j < geomL.nsurfz - 1; ++j) {

            int pid1L = geomL.surfp[i][j];
            int pid2L = geomL.surfp[i][j+1];
            
            int pid1R = geomR.surfp[i][j];
            int pid2R = geomR.surfp[i][j+1];

            if (pid1L < 0 || pid2L < 0 || pid1R < 0 || pid2R < 0) continue;

            // Z方向の幅 
            double yL_1 = stateL.disp[pid1L].uy ;
            double yL_2 = stateL.disp[pid2L].uy ;
            double yL_avg = 0.5 * (yL_1 + yL_2);

            // 右声帯のY座標 (現在の変位 + 速度×dt による予測)
            double yR_1 = stateR.disp[pid1R].uy ;
            double yR_2 = stateR.disp[pid2R].uy ;
            double yR_avg = 0.5 * (yR_1 + yR_2);

            double gap = yR_avg - yL_avg;

            // Z方向の幅 (現在の変位 + 速度×dt による予測)
            double zL_1 = stateL.disp[pid1L].uz ;
            double zL_2 = stateL.disp[pid2L].uz ;
            double ztmp = std::abs(zL_2 - zL_1);

            // 交差（めり込み）した場合は流路面積を0とする
            if (gap < 0.0) gap = 0.0;
            if (!std::isfinite(gap)) gap = 0.0;
            if (!std::isfinite(ztmp)) ztmp = 0.0;

            hi += gap * ztmp; 

            static int N = 0;

                // calcArea内のループの中で（最初の数ステップだけ）
            if ( N < 10) {
                std::cout << "Step: " << N
                        << " | yL_avg: " << yL_avg 
                        << " | yR_avg: " << yR_avg 
                        << " | gap: " << gap << std::endl;
                        N++;
            }
        }
        harea[i] = hi;
    }


    for (auto& anglePlane : degreeL) {
        for (auto& row : anglePlane) {
            std::fill(row.begin(), row.end(), 0.0);
        }
    }
    for (auto& anglePlane : degreeR) {
        for (auto& row : anglePlane) {
            std::fill(row.begin(), row.end(), 0.0);
        }
    }

    for (int i = 1; i < geomL.nxsup - 1; ++i) {
        for (int j = 1; j < geomL.nsurfz - 1; ++j) {
            
            int pid_leftL  = geomL.surfp[i-1][j];
            int pid_rightL = geomL.surfp[i+1][j];
            int pid_downL  = geomL.surfp[i][j-1];
            int pid_upL    = geomL.surfp[i][j+1];

            if (pid_leftL >= 0 && pid_rightL >= 0 && pid_downL >= 0 && pid_upL >= 0) {
                double dxL  = 0.5 * (stateL.disp[pid_rightL].ux - stateL.disp[pid_leftL].ux);
                double dy1L = 0.5 * (stateL.disp[pid_rightL].uy - stateL.disp[pid_leftL].uy);
                double dy2L = 0.5 * (stateL.disp[pid_upL].uy   - stateL.disp[pid_downL].uy);
                double dzL  = 0.5 * (stateL.disp[pid_upL].uz   - stateL.disp[pid_downL].uz);

                if (dxL != 0.0) degreeL[0][i][j] = std::atan(dy1L / dxL);
                if (dzL != 0.0) degreeL[1][i][j] = std::atan(dy2L / dzL);
            }


            int pid_leftR  = geomR.surfp[i-1][j];
            int pid_rightR = geomR.surfp[i+1][j];
            int pid_downR  = geomR.surfp[i][j-1];
            int pid_upR    = geomR.surfp[i][j+1];

            if (pid_leftR >= 0 && pid_rightR >= 0 && pid_downR >= 0 && pid_upR >= 0) {
                double dxR  = 0.5 * (stateR.disp[pid_rightR].ux - stateR.disp[pid_leftR].ux);
                double dy1R = 0.5 * (stateR.disp[pid_rightR].uy - stateR.disp[pid_leftR].uy);
                double dy2R = 0.5 * (stateR.disp[pid_upR].uy   - stateR.disp[pid_downR].uy);
                double dzR  = 0.5 * (stateR.disp[pid_upR].uz   - stateR.disp[pid_downR].uz);

                if (dxR != 0.0) degreeR[0][i][j] = std::atan(dy1R / dxR);
                if (dzR != 0.0) degreeR[1][i][j] = std::atan(dy2R / dzR);
            }
        }
    }
}

void ForceCalculator::calcDis() {

    for(int i = 0; i < geomL.nxsup; ++i){
        std::fill(contactForceL_ij[i].begin(), contactForceL_ij[i].end(), 0.0);
        std::fill(contactForceR_ij[i].begin(), contactForceR_ij[i].end(), 0.0);
    }


    contactFlag = false;

    // 固有角振動数 (非対称のため、左右の固有振動数の平均を剛性の基準とする)
    double omg1L = 2.0 * M_PI * modeDataL.frequencies[0];
    double omg1R = 2.0 * M_PI * modeDataR.frequencies[0];
    double omg2 = 0.5 * (omg1L * omg1L + omg1R * omg1R);

    max_force_diff = 0.0;

    for (int i = 1; i < nxsup; ++i) {
        for (int j = 1; j < geomL.nsurfz - 1; ++j) {
            
            // 前の反復で足した分を一度引いてキャンセルする
            fyL[i][j] -= fdisL[i][j]; 
            fyR[i][j] -= fdisR[i][j];

            double old_forceL = fdisL[i][j];
            double old_forceR = fdisR[i][j];

            fdisL[i][j] = 0.0;
            fdisR[i][j] = 0.0;

            int pidL = geomL.surfp[i][j];
            int pidR = geomR.surfp[i][j];

            if (pidL < 0 || pidR < 0) continue;

            // ループ内で更新された予測位置 
            double yL_curr = stateL.predictedDisp[pidL].uy; 
            double yR_curr = stateR.predictedDisp[pidR].uy; 

            if (yL_curr > yR_curr) {

                // --- 【追加】初回呼び出し時のみ実行される/* debug */コード ---
                static bool is_first_contact = true;
                if(is_first_contact){
                    std::cout << "\n=== [DEBUG] FIRST CONTACT DETECTED ===" << std::endl;
                    is_first_contact = false;
                }
                // 片側あたりのめり込み量
                double pen = 0.5 * (yL_curr - yR_curr) * 1e-3; // mm -> m

                // 2. 速度 (m/s) : (現在の予測位置 - 1ステップ前の位置) / dt
                double yL_prev = stateL.disp[pidL].uy;
                double velL = (yL_curr - yL_prev) / sp.dt * 1e-3; 

                double yR_prev = stateR.disp[pidR].uy;
                double velR = (yR_curr - yR_prev) / sp.dt * 1e-3; 

                double vel_rel = 0.5 * (velL - velR);

                // 3. バネ力（線形＋非線形）
                double non_linear_term = 1.0 + sp.kc2 * omg2 * pen * pen; 
                double f_spring = sp.kc1 * omg2 * pen * non_linear_term;
                double f_damp = sp.kc3 * vel_rel;

                // 4. 減衰力
                double contact_pressure = f_spring + f_damp;


                // 5. 合力 
                double areaL = geomL.sarea[i][j] * 1e-6;
                double areaR = geomR.sarea[i][j] * 1e-6;

                double f_totalL = -(f_spring + f_damp) * areaL;
                double f_totalR =  (f_spring + f_damp) * areaR;

                if (f_totalL > 0.0) f_totalL = 0.0;
                if (f_totalR < 0.0) f_totalR = 0.0;

                double diffL = std::abs(f_totalL - old_forceL);
                double diffR = std::abs(f_totalR - old_forceR);
                if (diffL > max_force_diff) max_force_diff = diffL;
                if (diffR > max_force_diff) max_force_diff = diffR;

                // 力を保存・適用
                fdisL[i][j] = f_totalL;
                contactForceL_ij[i][j] = f_totalL;
                fyL[i][j] += fdisL[i][j];

                fdisR[i][j] = f_totalR;
                contactForceR_ij[i][j] = f_totalR;
                fyR[i][j] += fdisR[i][j];

                contactFlag = true;
            }
        }
    }


}

void ForceCalculator::calcFlowStep(double t, double dt, double min_area) {

    static bool printed_constants = false;
    if (t > 0.0 && !printed_constants) {
        std::cout << "\n=== [DEBUG] calcFlowStep Constants at Step 1 ===" << std::scientific << std::setprecision(12) << std::endl;
        std::cout << "ps (Lung Press): " << sp.ps << std::endl;
        std::cout << "rho: " << rho << " | mu: " << mu << std::endl;
        std::cout << "xsup: " << geomL.xsup << " | lg: " << lg << std::endl;
        std::cout << "beta: " << beta << " | R2: " << R2 << std::endl;
        std::cout << "Cu: " << Cu << " | Lu: " << Lu << std::endl;
        std::cout << "La: " << La << " | Ca: " << Ca << std::endl;
        std::cout << "=================================================\n" << std::endl;
        printed_constants = true;
    }
    
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

double ForceCalculator::findMinHarea() {
    return *std::min_element(harea.begin(), harea.end());
}

int ForceCalculator::findNsep(double minH) {
    for (int i = 1; i < nxsup; i++) {
        if (std::fabs(harea[i] - minH) < 1e-8 || harea[i] <= 0.0) {
            return i+1;
        }
    }
    return geomR.nxsup;
}


void ForceCalculator::outputForceVectors(int step) const {
    std::ostringstream stepStr;
    stepStr << std::setw(4) << std::setfill('0') << step;

    std::ofstream fout("../output2/force_" + stepStr.str() + ".csv");
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
