#pragma once

#include <vector>
#include "ChannelSections.h"
#include "Geometry.h"
#include "SimulationParams.h"

// Time history and equations of the one-dimensional airway model.  It owns no
// surface mesh or structural state; geometry is supplied as ChannelSections.
class FlowModel {
public:
    void initialize(const SimulationParams& params, const Geometry& geometry,
                    int nSteps);
    void advance(double time, double dt, double minimumArea, double previousFlow,
                 const ChannelSections& sections);

    const std::vector<double>& pressure() const { return psurf_; }
    const std::vector<double>& area() const { return area_; }
    double flowRate() const { return currentUg_; }
    double glottalPressure() const { return currentPg_; }
    const std::vector<double>& upstreamPressure() const { return Pu_; }
    const std::vector<double>& downstreamPressure() const { return Pd_; }
    double outletPressure() const { return Pd_.empty() ? 0.0 : Pd_.back(); }

private:
    void calcFlowStep(double time, double dt, double minArea);
    int findSeparation(const ChannelSections& sections, double minArea) const;

    const SimulationParams* sp_ = nullptr;
    double rho_ = 0.0, mu_ = 0.0;
    double glottalInertanceGeometry_ = 0.0;
    double glottalViscousGeometry_ = 0.0;
    double alpha1_ = 0.0, alpha2_ = 0.0, beta_ = 0.0;
    double Lui_ = 0.0, Cui_ = 0.0, Lu_ = 0.0, Cu_ = 0.0, R2_ = 0.0;
    double La_ = 0.0, Ca_ = 0.0, Lr_ = 0.0, Rr_ = 0.0;
    bool hasVocalTract_ = false;
    double previousUg_ = 0.0, currentUg_ = 0.0, currentPg_ = 0.0;
    std::vector<double> Ug_, minAreaHistory_, Uu_, Pu_, Ud_, Pd_;
    std::vector<double> psurf_, area_;
};
