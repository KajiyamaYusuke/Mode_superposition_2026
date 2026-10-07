#include "SimulationParams.h"

#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: test_input_file_params PARAM_FILE\n";
        return 2;
    }

    SimulationParams params;
    std::string error;
    if (!params.loadFromFile(argv[1], error)) {
        std::cerr << error << "\n";
        return 1;
    }

    const bool matches =
        params.leftFrequencyFile == "M5_test/M5_freq_kawahara_mesh7_soft.txt"
        && params.rightFrequencyFile == "M5_test/M5_freq_kawahara_mesh7_soft.txt"
        && params.leftModeFile == "M5_test/M5_mode_kawahara_mesh7_soft.vtu"
        && params.rightModeFile == "M5_test/M5_mode_kawahara_mesh7_soft.vtu"
        && params.leftSurfaceNasFile == "M5_test/M5_surface_kawahara_mesh7_soft.nas"
        && params.rightSurfaceNasFile == "M5_test/M5_surface_kawahara_mesh7_soft.nas";
    if (!matches) {
        std::cerr << "fold-specific input file parameters were not parsed correctly\n";
        return 1;
    }
    if (params.apContactDiagnosticEnabled
        || params.apContactAnalysisRowsPerEnd != 3
        || params.apContactExcludedRowsPerEnd != 0
        || params.diagnosticLoadMode != "live"
        || params.diagnosticReferenceWindowStartSec != 0.15
        || params.diagnosticReferenceWindowEndSec != 0.19) {
        std::cerr << "AP contact diagnostic defaults were not parsed correctly\n";
        return 1;
    }
    return 0;
}
