#include <iostream>
#include <fstream>
#include <vector>

class TimeIntegrator {
public:

    
    void newmarkStep( double f, double q, double qdot, double qddot,
    double dt, double omg, double zeta,
    double beta, double gamma,
    double &qf, double &qfdot, double &qfddot);
};
