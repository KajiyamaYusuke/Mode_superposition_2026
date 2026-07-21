#include "FlowModel.h"

#include <algorithm>
#include <cmath>

void FlowModel::initialize(const SimulationParams& params, const Geometry& geometry,
                           int nSteps) {
    sp_ = &params;
    rho_ = params.rho; mu_ = params.mu; lg_ = geometry.zmax; xsup_ = geometry.xsup;
    Ug_.assign(nSteps, 0.0); minAreaHistory_.assign(nSteps, 0.0);
    Uu_.assign(params.N_sub + 1, 0.0);
    Pu_.assign(params.N_sub + 2, 0.0);
    Ud_.assign(params.N_vt, 0.0); Pd_.assign(params.N_vt, 0.0);

    psurf_.assign(geometry.nxsup, 0.0); area_.assign(geometry.nxsup, 0.0);

    const double Ain = M_PI*params.r_inlet*params.r_inlet;
    const double Asub = M_PI*params.r_sub*params.r_sub;
    const double Avt = M_PI*params.r_vt*params.r_vt;

    Lui_ = rho_*params.L_inlet/(2.0*Ain); Cui_ = params.L_inlet*Ain/(rho_*params.c_sound*params.c_sound);

    const double dxSub = params.L_sub/std::max(1, params.N_sub);

    Lu_ = rho_*dxSub/(2.0*Asub); Cu_ = dxSub*Asub/(rho_*params.c_sound*params.c_sound);

    alpha1_ = -2.5e-5*params.ps+0.185; alpha2_ = 1.6e-3*params.ps+0.6;
    beta_ = 1.125e-4*params.ps+0.1375;
    R2_ = alpha1_/(Asub*Asub)*std::sqrt(rho_*mu_*params.c_sound);
    
    hasVocalTract_ = params.L_vt > 1e-6 && Avt > 1e-12;
    if (hasVocalTract_) {
        const double dxVt = params.L_vt/std::max(1, params.N_vt);
        La_ = rho_*dxVt/Avt; Ca_ = dxVt*Avt/(rho_*params.c_sound*params.c_sound);
        Lr_ = rho_*1.1*std::sqrt(Avt/M_PI)/Avt;
        Rr_ = alpha2_*rho_*params.c_sound/(9.0*M_PI*M_PI*Avt);
    } else { La_ = 1e-1; Ca_ = 1e30; Lr_ = Rr_ = 0.0; }
}

void FlowModel::advance(double time, double dt, double minimumArea, double previousFlow) {
    previousUg_ = previousFlow;
    calcFlowStep(time, dt, minimumArea);
}

void FlowModel::calcFlowStep(double time, double dt, double minimumArea) {
    const int ng = sp_->N_sub, nv = sp_->N_vt;
    const double ramp = time < 0.05 ? 0.5*(1.0-std::cos(M_PI*time/0.05)) : 1.0;
    const double lungPressure = sp_->ps*ramp;
    for (int j=0; j<ng; ++j) Pu_[j] += Uu_[j]-Uu_[j+1];
    Pu_[ng] += Uu_[ng]-previousUg_;
    Pu_[ng+1] += previousUg_-Ud_[0];
    Uu_[0] -= dt/Lui_*(dt/Cui_*Pu_[0]-lungPressure);
    Uu_[1] -= dt/(Lui_+Lu_)*(dt/Cu_*Pu_[1]-dt/Cui_*Pu_[0]+R2_*Uu_[1]);
    for (int j=2; j<ng+1; ++j)
        Uu_[j] -= dt/(2.0*Lu_)*(dt/Cu_*Pu_[j]-dt/Cu_*Pu_[j-1]+R2_*Uu_[j]);
    if (minimumArea > 1e-8) {
        const double lis = xsup_*1e-3;
        const double Lg1 = rho_*0.5*lis/minimumArea;
        const double Rk1 = beta_*rho_/(minimumArea*minimumArea);
        const double Rv1 = 12.0*mu_*lis*(lg_*1e-3)*(lg_*1e-3)/std::pow(minimumArea,3.0);
        double guess = currentUg_;
        for (int k=0;k<100;++k) {
            const double f = Rk1*std::abs(guess)*guess + Rv1*guess
                +(Lg1+La_+Lu_)*(guess-previousUg_)/dt + dt/Ca_*Pu_[ng+1]-dt/Cu_*Pu_[ng];
            const double df = 2.0*Rk1*std::abs(guess)+Rv1+(Lg1+La_+Lu_)/dt;
            if (std::abs(f)<1e-9) break;
            guess -= f/df;
        }
        currentUg_ = guess;
    } else currentUg_ = 0.0;
    currentPg_ = dt/Cu_*Pu_[ng];
    if (!hasVocalTract_) {
        for (int i=0;i<nv;++i) { Pd_[i]=0.0; Ud_[i]=currentUg_; }
    } else {
        Pd_[0] += dt/Ca_*(currentUg_-Ud_[0]);
        for(int i=1;i<nv;++i) Pd_[i] += dt/Ca_*(Ud_[i-1]-Ud_[i]);
        for(int i=0;i<nv-1;++i) Ud_[i] += dt/La_*(Pd_[i]-Pd_[i+1]);
    }
    const double zrad = (La_+Lr_)/dt+Rr_;
    Ud_[nv-1] = (Pd_[nv-1]+(La_+Lr_)/dt*Ud_[nv-1])/zrad;
}
