#pragma once

#include <memory>

class Simulation;

// Read-only observer. Compiled only with ENABLE_KANAYA_DIAGNOSTICS, and
// inactive unless SIMULATION_KANAYA_DIAGNOSTICS=1. Never owns solver state.
class KanayaDiagnostics {
public:
    explicit KanayaDiagnostics(Simulation& simulation);
    ~KanayaDiagnostics();
    void beginStep(int step, double time);
    void afterFluid();
    void adoptedContact();
    void endStep(int iterations, int stopReason);
    void finish();

    static void node(bool left, int i, int j, double pressure,
                     double ax, double ay, double az,
                     double fx, double fy, double fz);
    static void beginContact(bool diagnosticOnly, int iteration);
    static void contact(int j, double x, double z, double penetrationM,
                        double fxL, double fyL, double fxR, double fyR,
                        bool excluded);
    static void coefficients(double lg, double rv, double rk, bool active);
    static double openingForce(double gapSign, bool left, double fy);
    static double scaledResidual(double difference, double scale,
                                 double absoluteFloor = 1.0e-12);
    static bool includedInterval(int step, int nSteps);

private:
    struct Data;
    std::unique_ptr<Data> data_;
    static Data* active_;
};
