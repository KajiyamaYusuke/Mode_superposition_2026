#include "Simulation.h"
#include <iostream>
#include <iomanip>
#include <iostream>

int main() {
    Simulation sim;

    sim.initialize();
    sim.params.iforce = 0;

    sim.run();


    return 0;
}



