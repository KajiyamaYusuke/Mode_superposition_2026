#include "SimulationParams.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void expectNear(double actual, double expected, const char* name) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-8) {
        throw std::runtime_error(std::string(name) + ": unexpected pressure");
    }
}
void expectValid(const SimulationParams& p) {
    std::string error;
    if (!p.validate(error)) throw std::runtime_error(error);
}
void expectInvalid(const SimulationParams& p) {
    std::string error;
    if (p.validate(error)) throw std::runtime_error("expected invalid pressure settings");
}
}
int main() {
    try {
        SimulationParams legacy;
        expectValid(legacy);
        expectNear(legacy.pressureCommandAt(0.0), 0.0, "legacy start");
        expectNear(legacy.pressureCommandAt(0.025), legacy.ps * 0.5, "legacy midpoint");
        expectNear(legacy.pressureCommandAt(0.05), legacy.ps, "legacy end");

        SimulationParams rise;
        rise.pressureProfile = "linear_sweep";
        rise.pressureInitialPa = 0.0;
        rise.pressureTargetPa = 2000.0;
        rise.pressureRiseTimeSec = 1.0;
        rise.pressureInitialPaExplicit = true;
        rise.pressureTargetPaExplicit = true;
        rise.pressureRiseTimeSecExplicit = true;
        expectValid(rise);
        expectNear(rise.pressureCommandAt(0.0), 0.0, "rise start");
        expectNear(rise.pressureCommandAt(0.5), 1000.0, "rise midpoint");
        expectNear(rise.pressureCommandAt(1.0), 2000.0, "rise end");
        expectNear(rise.pressureCommandAt(1.4), 2000.0, "rise hold");

        SimulationParams fall = rise;
        fall.pressureInitialPa = 2000.0;
        fall.pressureTargetPa = 0.0;
        expectValid(fall);
        expectNear(fall.pressureCommandAt(0.0), 2000.0, "fall start");
        expectNear(fall.pressureCommandAt(0.5), 1000.0, "fall midpoint");
        expectNear(fall.pressureCommandAt(1.0), 0.0, "fall end");

        SimulationParams cycle = rise;
        cycle.pressureProfile = "linear_cycle";
        cycle.pressureFinalPa = 0.0;
        cycle.pressureFinalPaExplicit = true;
        cycle.pressureHoldTimeSec = 0.2;
        cycle.pressureHoldTimeSecExplicit = true;
        cycle.pressureFallTimeSec = 1.0;
        cycle.pressureFallTimeSecExplicit = true;
        expectValid(cycle);
        expectNear(cycle.pressureCommandAt(0.5), 1000.0, "cycle rise");
        expectNear(cycle.pressureCommandAt(1.1), 2000.0, "cycle peak");
        expectNear(cycle.pressureCommandAt(1.7), 1000.0, "cycle fall");
        expectNear(cycle.pressureCommandAt(2.2), 0.0, "cycle final");

        SimulationParams invalid = rise;
        invalid.pressureRampTimeSecExplicit = true;
        expectInvalid(invalid);
        invalid = rise;
        invalid.pressureRiseTimeSec = 0.0;
        expectInvalid(invalid);
        invalid = rise;
        invalid.pressureTargetPa = 8000.0;
        expectInvalid(invalid);
        invalid = rise;
        invalid.iforce = 1;
        expectInvalid(invalid);

        std::cout << "pressure profile tests passed\n";
    } catch (const std::exception& e) {
        std::cerr << "pressure profile test failed: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
