#include "TimeIntegrator.h"
#include <vector>


void TimeIntegrator::newmarkStep(
    double f, double q, double qdot, double qddot,
    double dt, double omg, double zeta,
    double beta, double gamma,
    double &qf, double &qfdot, double &qfddot)
{
    double c = 2.0 * zeta * omg;   // damping
    double k = omg * omg;          // stiffness

    // Predictor
    double q_pred    = q + dt * qdot + dt * dt * (0.5 - beta) * qddot;
    double qdot_pred = qdot + dt * (1.0 - gamma) * qddot;

    // Effective stiffness (denominator)
    double keff = 1.0 + gamma * dt * c + beta * dt * dt * k;

    // Effective force
    double feff = f
        - k * q_pred
        - c * qdot_pred;

    // Solve for new acceleration
    double qdd_new = feff / keff;

    // Corrector
    qf     = q_pred    + beta  * dt * dt * qdd_new;
    qfdot  = qdot_pred + gamma * dt      * qdd_new;
    qfddot = qdd_new;
}

