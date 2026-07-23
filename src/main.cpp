#include "Simulation.h"
#include <iostream>
#include <filesystem>

int main(int argc, char* argv[]) {
    Simulation sim;

    const fs::path parameterFile = argc >= 2 ? argv[1] : "../input/param.txt";
    sim.initialize(parameterFile);

    sim.run();


    return 0;
}


