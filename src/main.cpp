#include "Simulation.h"
#include <iostream>
#include <filesystem>

int main(int argc, char* argv[]) {
    try {
        Simulation sim;
        const fs::path parameterFile = argc >= 2 ? argv[1] : "../input/param.txt";
        sim.initialize(parameterFile);
        sim.run();
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "[Simulation] ERROR: " << exception.what() << "\n";
        return 1;
    }
}

