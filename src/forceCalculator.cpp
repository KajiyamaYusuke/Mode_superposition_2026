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



ForceCalculator::ForceCalculator(const Geometry& geom_, const ModeData& md_, State& st_, const SimulationParams& sp_)
    : geom(geom_), modeData(md_), state(st_), sp(sp_) {}




void ForceCalculator::initialize() {
    int nPoints = geom.nPoints;
    int nsurfl  = geom.nsurfl;
    int nsurfz  = geom.nsurfz;
    int nModes  = modeData.nModes;

    // 節点ごとの外力
    fx.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fy.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fz.assign(nsurfl, std::vector<double>(nsurfz, 0.0));

    // モード力
    fi.assign(nModes, 0.0);

    // 内部バッファ
    nxsup = geom.nxsup;  // 断面数 (x方向分割)
    fdis.assign(nxsup, std::vector<double>(geom.nsurfz - 1, 0.0));

    psurf.assign(nxsup, 0.0);
    Ug.assign(state.nSteps, 0.0);
    minHarea.assign(state.nSteps, 0.0);

    Uu.assign(Nsecg + 1, 0.0); // +2 for boundaries
    Pu.assign(Nsecg + 2, 0.0);
    Ud.assign(Nsecp , 0.0);
    Pd.assign(Nsecp , 0.0);

    
    currentUg = 0.0;    

    rho = sp.rho ;
    mu  = sp.mu ;
    lg  = geom.zmax ; 

    c_sound = sp.c_sound;

    double r_chamber = 16.0 * 1e-2; // 半径 15cm -> m
    double A_inlet = M_PI * r_chamber * r_chamber;
    double L_inlet = 50.0  * 1e-2; // cm -> m
    
    // 2. Subglottal Tract (声門下)
    double A_sub   = M_PI * std::pow((1.6 * 1e-2) / 2.0, 2.0);
    double L_sub   = 46.0  * 1e-2;
    int    N_sub   = Nsecg; // param.txt の section数

    // 3. Vocal Tract (声道)
    double r_tract = 0.8 * 1e-2; // 半径 15cm -> m
    double A_vt  = M_PI * r_tract * r_tract;
    double L_vt    = 17.0 * 1e-2;
       int N_vt    = 10; // param.txt の section数 (Nsecpで使用)

    // --- インピーダンス計算 (L = rho*l/A, C = V / (rho*c^2) = l*A / (rho*c^2)) ---
    
    // Inlet parameters (Lui, Cui)
    // Fortranでは Lui = rho * L_inlet / A_inlet など
    Lui = rho * L_inlet / (2*A_inlet);
    Cui = L_inlet * A_inlet / (rho * c_sound * c_sound);

    // Subglottal parameters (Lu, Cu) - 1セクションあたり
    double dx_sub = L_sub / N_sub;
    Lu = rho * dx_sub / A_sub;
    Cu = dx_sub * A_sub / (rho * c_sound * c_sound);

    alpha1 = -2.5e-5*sp.ps+0.185;
    alpha2 = 1.6e-3*sp.ps+0.6;
    beta = 1.125e-4 * sp.ps + 0.1375;

    R2 = alpha1/(A_vt*A_vt) * std::sqrt(rho*mu*c_sound);
    if (L_vt > 1e-6) {
        double dx_vt = L_vt / std::max(1, Nsecp); 
        La = rho * dx_vt / A_vt;
        Ca = dx_vt * A_vt / (rho * c_sound * c_sound);
        
        Lr = rho * 1.1 * sqrt(A_vt/M_PI) /A_vt;
        Rr = alpha2*rho*c_sound/(9*M_PI*M_PI*A_vt);
    } else {
        // 声道がない場合、計算に使わないがゼロ除算回避のため安全な値を入れておく
        // または calcFlowStep で分岐する
        La = 1e-1; // ダミー
        Ca = 1.0e30; // 非常に大きくすることで圧力変動をゼロにする(大気開放)
        Lr = 0.0;
        Rr = 0.0;
    }
    

    std::cout << "[ForceCalculator] initialized: "
              << "nPoints=" << nPoints
              << ", nModes=" << nModes
              << ", nxsup=" << nxsup 
              << ", L_sub=" << L_sub
              << ", L_vt=" << L_sub
              << std::endl;
}

void ForceCalculator::calcForce(double t, int n) {
    int nsurfl = geom.nsurfl;
    int nsurfz = geom.nsurfz;
    int nxsup  = geom.nxsup;

    // まず全てゼロクリア
    fx.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fy.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fz.assign(nsurfl, std::vector<double>(nsurfz, 0.0));


    for (int i = 0; i < nxsup; i++) {
        for (int j = 0; j < nsurfz - 1; j++) {
            fdis[i][j] = 0.0;
        }
    }


    if (sp.iforce == 1) {
        // ==== sin波加振 ====

        minHarea[n] = *std::min_element(state.harea.begin(), state.harea.end());
        
        
        for (int i = 1; i < nxsup-1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {

                fx[i][j] = sp.famp * std::sin(2.0* M_PI *sp.forcef* t);
                //fy[i][j] = sp.famp * std::sin(2.0* M_PI *sp.forcef* t);
                //fz[i][j] = sp.famp * std::sin(2.0* M_PI *sp.forcef* t);
            }
        }
    } else if (sp.iforce == 0) {

        // ==== 1D flow model ====
        double minA = findMinHarea();
        minHarea[n] = minA;

        previousUg = (n > 0 && n-1 < (int)Ug.size()) ? Ug[n-1] : 0.0;
        calcFlowStep(t, sp.dt, minA * 1e-6);
        // separation point
        int nsep = findNsep(minA) ;

        Ug[n] = currentUg;

        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = currentPg;

        if (minA > 1e-6 ) {
            for (int i = 1; i < geom.nxsup; i++) {
                double dx = std::abs(geom.points[geom.surfp[i][ int(nsurfz/2)]].x - geom.points[geom.surfp[i-1][int(nsurfz/2)]].x);
                double h  = (state.harea[i] + state.harea[i-1]) / (2.0 * lg);
                double h_prev = std::max(state.harea[i-1], 1e-6);
                double h_curr = std::max(state.harea[i], 1e-6);

                double Ugm = currentUg*1e6;

                double bernoulli_term = 0.5 * rho * Ugm * Ugm * (1.0/(h_prev*h_prev) - 1.0/(h_curr*h_curr));
                double viscous_term   = 12.0 * mu * dx * Ugm / (lg * pow(h, 3.0)) *1e3;

                psurf[i] = psurf[i-1] + bernoulli_term - viscous_term;

                if(psurf[i] > psurf[i-1]){break;}

            }


        } else {
            for (int i = 1; i < nsep-1; ++i) psurf[i] = currentPg;
            for (int i = nsep-1; i < nxsup; ++i) psurf[i] = Pd[0];
        }


        // 力 fx, fy, fz
        for (int i = 1; i < nsep-1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {
                int pid_ip1 = geom.surfp[i+1][j];
                int pid_im1 = geom.surfp[i-1][j];
                int pid_jp1 = geom.surfp[i][j+1];
                int pid_jm1 = geom.surfp[i][j-1];

                double dx = 0.5 * (state.disp[pid_ip1].ux - state.disp[pid_im1].ux);
                double dy = 0.5 * (state.disp[pid_ip1].uy - state.disp[pid_im1].uy);
                double ds = std::sqrt(dx*dx + dy*dy);
                double dz = 0.5 * (state.disp[pid_jp1].uz - state.disp[pid_jm1].uz);

                fx[i][j] = psurf[i] * ds * dz * 1.0e-6 * std::cos(state.degree[1][i][j]) * std::sin(state.degree[0][i][j]);
                fy[i][j] = -psurf[i] * ds * dz * 1.0e-6 * std::cos(state.degree[1][i][j]) * std::cos(state.degree[0][i][j]);
                fz[i][j] = psurf[i] * ds * dz * 1.0e-6 * std::sin(state.degree[1][i][j]);

            }
        }
        
    }
    

}

void ForceCalculator::f2mode() {



    for (int imode = 0; imode < modeData.nModes; imode++) {
        fi[imode] = 0.0;
        for (int i = 0; i < geom.nsurfl; i++) {
            for (int j = 0; j < geom.nsurfz; j++) {
                int pid = geom.surfp[i][j];
                if (pid < 0) continue; // 節点が存在しない場合スキップ

                fi[imode] += fx[i][j] * modeData.modes[imode][pid].ux
                           + fy[i][j] * modeData.modes[imode][pid].uy
                           + fz[i][j] * modeData.modes[imode][pid].uz;

            }
        }
    }

     //std::cout<<"|disp[1]= "<<state.disp[11].ux<<std::endl;
}

void ForceCalculator::contactForce() {

    contactFlag = false;
    double omg1 = 2.0 * M_PI * modeData.frequencies[0]; // 1次固有振動数
    double omg2 = omg1 * omg1;

    for (int i = 0; i < geom.nxsup; ++i) {           // nxsup は計算範囲
        for (int j = 1; j < geom.nsurfz - 1; ++j) {  // 2..nsurfz-1 (0-index)
            int node = geom.surfp[i][j];

            double y     = state.disp[node].uy;        // 現在変位
            double ydot  = state.vel[node].uy;         // 現在速度
            double ymid  = geom.ymid[j];
            double yhat  = y + sp.dt * ydot;           // 予測位置


            // 接触状態を判定
            bool contact_now    = (y > ymid);          // 現時点で接触
            bool contact_future = (yhat > ymid);       // 次ステップで接触

            if (!contact_now) {
                continue;
            }

            double pen = (ymid - y) * 1e-3; 

            double f_contact = sp.kc1 * omg2 * pen * (1.0 + sp.kc2 * omg2 * pen * pen);

            double f_damp =  sp.kc3 * pen * ydot;

            double f_total = (f_contact + f_damp) * geom.sarea[i][j] * 1e-6;

            if (f_total > 0.0) { f_total = 0.0; }           

            fy[i][j] += f_total;
            contactFlag = true;
        }
    } 

}

void ForceCalculator::calcDis() {

    contactFlag = false;

    // 単位質量の仮定
    double mass = sp.mass;

    for (int i = 1; i < nxsup; ++i) {            // 2..nxsup (0-indexなので1スタート)
        for (int j = 1; j < geom.nsurfz - 1; ++j) { 
            fdis[i][j] = 0.0;

            int pid = geom.surfp[i][j];
            if (pid < 0) continue;

            double v_now = state.disp[pid].uy;
            double v_next = state.predictedDisp[pid].ufy;           // Runge-Kuttaで予測した変位
            

            if (v_now <= geom.ymid[j] && v_next > geom.ymid[j]) {
                double vc = (v_next - v_now) / sp.dt; // mm/s -> m/s は必要なら換算
                vc *= 1e-3;

                fdis[i][j] = -mass * vc / (sp.dt * geom.nPoints);
                fy[i][j] += fdis[i][j];

                contactFlag = true;
            }
        }
    }
}

void ForceCalculator::calcFlowStep(double t, double dt, double min_area) {
    
    // --- 1. 声門下 (Subglottal) の更新 ---
    
    double rampDuration = 0.075; // 50msかけて立ち上げる
    double rampFactor = 1.0;
    
    if (t < rampDuration) {
        // Cosine Ramp (滑らか)
        rampFactor = 0.5 * (1.0 - std::cos(M_PI * t / rampDuration));
    }

    // --- ランプ適用 ---
    // sp.ps (固定パラメータ) に rampFactor をかけて「現在の肺圧」を作る
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


    // 流量の更新 (運動量保存: dU/dt = (1/L) * (Pin - Pout - R*U))
    // Uu[1]: Inlet -> 1st Section
    // Fortran: Uu(1)=Uu(1)-dt/Lui*(dt/Cui*Pu(1)-Ps) 
    // これは「P(1) - Ps」の形。
    // C++:
    Uu[0]  -= dt / Lui * ( (dt / Cui * Pu[0]) - currentLungPressure );

    // Fortran: Uu(2)=Uu(2)-dt/(Lui+Lu)*(dt/Cu*Pu(2)-dt/Cui*Pu(1)+R2*Uu(2))
    Uu[1] -= (dt / (Lui + Lu)) * ( dt / Cu * Pu[1] - dt / Cui * Pu[0] + R2 * Uu[1] );

    for (int j = 2; j < Nsecg + 1; ++j) {
        Uu[j] -= dt / (2.0 * Lu) * ( dt / Cu * Pu[j] - dt / Cu * Pu[j-1] + R2 * Uu[j] );
    }


    // --- 2. 声門部 (Glottal Flow) の更新 ---
    if (min_area > 1e-8) {
        double min_area_m2 = min_area;
        double lis = geom.xsup * 1e-3; // 仮定値 (Fortran側でd1+d2に相当するか要確認)
        double lg_m = lg * 1e-3;

        double Lg1 = rho *  0.5 * lis / min_area_m2;
        double Rk1 = beta * rho / ( min_area_m2 * min_area_m2); // Bernoulli (係数調整)
        // Fortranでは beta*rho... とある。betaが1.0以上なら損失
        double Rv1 = 12.0 * mu * lis * lg_m * lg_m / pow(min_area_m2, 3.0);

        // 駆動圧: 声門直下(Pu[last]) - 声道入口(Pd[0])
        double Ug_old = previousUg; // これを適切に保持しておく
        double Ug_guess = currentUg; // or previous guess

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

    // currentPg (声門下圧として外力計算に使う値)
    currentPg = dt / Cu * Pu[Nsecg];


    // --- 3. 声道 (Vocal Tract) の更新 ---
    // 圧力更新 Pd[0]...

    if(Lr < 1e-8){
        for (int i = 0; i < Nsecp; i++) {
            Pd[i] = 0.0;     // 圧力ゼロ（大気）
            Ud[i] = currentUg; // 流量はすべて Ug と同じ
        }
    }else{
        Pd[0] += (dt / Ca) * (currentUg - Ud[0]);
        for(int i=1; i<Nsecp; ++i) {
            Pd[i] += (dt / Ca) * (Ud[i-1] - Ud[i]);
        }

        // 流量更新 Ud[0]...
        for(int i=0; i<Nsecp-1; ++i) {
            Ud[i] += (dt / La) * (Pd[i] - Pd[i+1]); // 抵抗Raがあれば追加
        }
    }
    
    // 放射端 (Radiation)
    // Lr * dUd/dt + Rr * Ud = Pd[last] - P_atm(0)
    // 離散化: (Lr/dt + Rr) * Ud_new = Pd[last] + (Lr/dt)*Ud_old
    double Z_rad = (La+Lr)/dt + Rr;
    Ud[Nsecp-1] = (Pd[Nsecp-1] + ((La+Lr)/dt)*Ud[Nsecp-1]) / Z_rad;
}

double ForceCalculator::findMinHarea() {
    return *std::min_element(state.harea.begin(), state.harea.end());
}

int ForceCalculator::findNsep(double minH) {
    for (int i = 1; i < geom.nxsup; i++) {
        if (std::fabs(state.harea[i] - minH) < 1e-8 || state.harea[i] <= 0.0) {
            return i+1;
        }
    }
    return geom.nxsup;
}


void ForceCalculator::outputForceVectors(int step) const {
    std::ostringstream stepStr;
    stepStr << std::setw(4) << std::setfill('0') << step;

    std::ofstream fout("../output2/force_" + stepStr.str() + ".csv");
    fout << "x,y,z,Fx,Fy,Fz\n";  // CSVヘッダー

    for (int i = 1; i < geom.nxsup - 1; ++i) {
        for (int j = 1; j < geom.nsurfz - 1; ++j) {
            int pid = geom.surfp[i][j];
            if (pid < 0) continue;

            const auto &p = geom.points[pid];
            fout << p.x << "," << p.y << "," << p.z << ","
                 << fx[i][j] << "," << fy[i][j] << "," << fz[i][j] << "\n";
        }
    }

    //std::cout << "[Output] force vectors written for step " << step << std::endl;
}
