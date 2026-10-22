#include "FlowModel.h"
#ifdef ENABLE_KANAYA_DIAGNOSTICS
#include "KanayaDiagnostics.h"
#endif
#include "FluidConstants.h"
#include "SeparationSelection.h"
#include "LegacyGlottalMethods.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstdlib>

void FlowModel::initialize(const SimulationParams& params, const Geometry& geometry,
                           int nSteps) {
    sp_ = &params;
    const char* trace = std::getenv("SIMULATION_ACOUSTIC_TRACE");
    traceAcoustics_ = params.upstreamAcousticProfile == "legacy_main"
        || (trace && std::string(trace) == "1");
    rho_ = params.rho; mu_ = params.mu;
    if (params.glottalMethods.any()) {
        double minX=geometry.points.front().x, maxX=minX;
        double minZ=geometry.points.front().z, maxZ=minZ;
        for (const auto& p:geometry.points) {
            minX=std::min(minX,p.x); maxX=std::max(maxX,p.x);
            minZ=std::min(minZ,p.z); maxZ=std::max(maxZ,p.z);
        }
        // The user's successful C++ uses absolute xsup/zmax, whereas public
        // Fortran uses the x range for lis. They agree on the supplied mesh.
        legacyLisM_=maxX*1e-3;
        legacyLgM_=geometry.zmax*1e-3;
    }
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

    if (params.upstreamAcousticProfile == "legacy_main")
        Lu_ = rho_*dxSub/Asub;
    else
        Lu_ = rho_*dxSub/(2.0*Asub);
    Cu_ = dxSub*Asub/(rho_*params.c_sound*params.c_sound);

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

void FlowModel::advance(double time, double dt, double minimumArea, double previousFlow,
                        const ChannelSections& sections) {
    const auto selection = sp_->glottalMethods.any()
        ? LegacyGlottal::select(sections,sp_->glottalMethods,sp_->flowSeparationAreaRatio)
        : selectSeparationByAreaRatio(sections.area, sections.valid, sp_->flowSeparationAreaRatio, kAreaClosedMm2);
    auto fail = [time](int i, const std::string& reason) {
        throw std::runtime_error("flow time_s=" + std::to_string(time)
            + " section=" + std::to_string(i) + ": " + reason);
    };
    if (selection.sepIndex < 0) fail(selection.minIndex, selection.reason);
    if (sections.x.size() != sections.area.size()
        || sections.gapCubed.size() != sections.area.size())
        fail(-1, "inconsistent section array sizes");
    if (minimumArea != selection.minAreaMm2 * 1.0e-6)
        fail(selection.minIndex, "inconsistent minimum area");
    if (sp_->glottalMethods.flowCoefficients) {
        glottalInertanceGeometry_=glottalViscousGeometry_=0.;
        if (minimumArea>kAreaClosedM2) {
            const auto old=LegacyGlottal::coefficients(rho_,mu_,beta_,legacyLisM_,legacyLgM_,minimumArea);
            glottalInertanceGeometry_=old.Lg/rho_;
            glottalViscousGeometry_=old.Rv/(12.*mu_);
        }
        previousUg_=previousFlow;
        calcFlowStep(time,dt,minimumArea);
        if (traceAcoustics_) recordAcousticStep(dt,minimumArea);
        return;
    }
    if (sp_->glottalMethods.area && minimumArea<=kAreaClosedM2) {
        glottalInertanceGeometry_=glottalViscousGeometry_=0.;
        previousUg_=previousFlow;
        calcFlowStep(time,dt,minimumArea);
        if (traceAcoustics_) recordAcousticStep(dt,minimumArea);
        return;
    }
    // Distributed glottal inertance and viscosity are calculated from the
    // same selection as wall pressure. Both Lg and Rv stop at separation.
    // Only adjacent valid sections contribute; never bridge an invalid cut.
    double invAreaIntegral = 0.0; // [1/m]
    double invGapCubedIntegral = 0.0; // [1/m^3]
    for (int i = 0; i <= selection.sepIndex; ++i) {
        if (!sections.valid[i]) continue;
        if (!std::isfinite(sections.x[i]) || !std::isfinite(sections.area[i])
            || sections.area[i] < 0.0 || !std::isfinite(sections.gapCubed[i])
            || sections.gapCubed[i] < 0.0) fail(i, "invalid geometric integral input");
        if (sections.area[i] > kAreaClosedMm2 && !(sections.gapCubed[i] > 0.0))
            fail(i, "open section has nonpositive G3");
        if (i == 0 || !sections.valid[i - 1]) continue;
        const double dxMm = sections.x[i] - sections.x[i - 1];
        const double areaMm2 = 0.5 * (sections.area[i] + sections.area[i - 1]);
        const double g3Mm4 = 0.5 * (sections.gapCubed[i] + sections.gapCubed[i - 1]);
        if (!(dxMm > 0.0) || !std::isfinite(dxMm)
            || !std::isfinite(areaMm2) || !std::isfinite(g3Mm4))
            fail(i, "nonmonotonic x or nonfinite section average");
        if (areaMm2 > kAreaClosedMm2) invAreaIntegral += dxMm * 1.0e3 / areaMm2;
        if (g3Mm4 > 1.0e-12)
            invGapCubedIntegral += dxMm * 1.0e9 / g3Mm4;
    }
    if (!std::isfinite(invAreaIntegral) || !std::isfinite(invGapCubedIntegral))
        fail(selection.sepIndex, "nonfinite Lg/Rv geometry");
    previousUg_ = previousFlow;
    // Store the equivalent geometric factors in the existing members.  They
    // are not independent acoustic state and are rebuilt every time step.
    glottalInertanceGeometry_ = 0.5 * invAreaIntegral;
    glottalViscousGeometry_ = invGapCubedIntegral;
    calcFlowStep(time, dt, minimumArea);
    if (traceAcoustics_) recordAcousticStep(dt,minimumArea);
}

void FlowModel::calcFlowStep(double time, double dt, double minimumArea) {
#ifdef ENABLE_KANAYA_DIAGNOSTICS
    KanayaDiagnostics::coefficients(0.0, 0.0, 0.0, false);
#endif
    const int ng = sp_->N_sub, nv = sp_->N_vt;
    // The checked-in solver historically used 0.05 s although the model
    // documentation stated 0.10 s.  Preserve omitted-option runs exactly;
    // new diagnostic inputs should set pressureRampTimeSec explicitly.
    const double rampTime = sp_->pressureRampTimeSecExplicit
        ? sp_->pressureRampTimeSec : 0.05;
    rampRatio_ = time < rampTime
        ? 0.5 * (1.0 - std::cos(M_PI * time / rampTime)) : 1.0;
    lungPressure_ = sp_->ps * rampRatio_;
    for (int j=0; j<ng; ++j) Pu_[j] += Uu_[j]-Uu_[j+1];
    Pu_[ng] += Uu_[ng]-previousUg_;
    Pu_[ng+1] += previousUg_-Ud_[0];
    Uu_[0] -= dt/Lui_*(dt/Cui_*Pu_[0]-lungPressure_);
    Uu_[1] -= dt/(Lui_+Lu_)*(dt/Cu_*Pu_[1]-dt/Cui_*Pu_[0]+R2_*Uu_[1]);
    for (int j=2; j<ng+1; ++j)
        Uu_[j] -= dt/(2.0*Lu_)*(dt/Cu_*Pu_[j]-dt/Cu_*Pu_[j-1]+R2_*Uu_[j]);
    if (minimumArea > kAreaClosedM2) {
        const double Lg1 = rho_ * glottalInertanceGeometry_;
        const double Rk1 = beta_*rho_/(minimumArea*minimumArea);
        const double Rv1 = 12.0 * mu_ * glottalViscousGeometry_;
#ifdef ENABLE_KANAYA_DIAGNOSTICS
        KanayaDiagnostics::coefficients(Lg1, Rv1, Rk1, true);
#endif
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

// Keep the production NR kernel unchanged. Observing its intermediate
// expressions inside the loop can change compiler contraction/rounding.
// This separate post-step observer never writes an acoustic state.
[[gnu::noinline]] void FlowModel::recordAcousticStep(double dt, double minimumArea) {
    newtonTrace_ = NewtonTrace{};
    newtonTrace_.previousFlow = previousUg_;
    if (minimumArea <= kAreaClosedM2) return;
    newtonTrace_.active = true;
    newtonTrace_.iterations = -1; // Deliberately not instrumented inside NR.
    const double lg=rho_*glottalInertanceGeometry_;
    const double rk=beta_*rho_/(minimumArea*minimumArea);
    const double rv=12.0*mu_*glottalViscousGeometry_;
    newtonTrace_.terms = {rk*std::abs(currentUg_)*currentUg_,rv*currentUg_,
        (lg+La_+Lu_)*(currentUg_-previousUg_)/dt,
        dt/Ca_*Pu_[sp_->N_sub+1],-dt/Cu_*Pu_[sp_->N_sub]};
    newtonTrace_.residual = newtonTrace_.terms[0]+newtonTrace_.terms[1]
        +newtonTrace_.terms[2]+newtonTrace_.terms[3]+newtonTrace_.terms[4];
    // The original solver target is 1e-9 Pa. Independent recomputation has
    // a 1e-8 Pa roundoff allowance, identical to the previous analysis gate.
    newtonTrace_.converged = std::isfinite(newtonTrace_.residual)
        && std::abs(newtonTrace_.residual)<1e-8;
    if (sp_->upstreamAcousticProfile=="legacy_main" && !newtonTrace_.converged)
        throw std::runtime_error("S3 post-step NR residual failed: " + std::to_string(newtonTrace_.residual));
}
