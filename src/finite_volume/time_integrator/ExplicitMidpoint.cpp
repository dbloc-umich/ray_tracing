#include "ExplicitMidpoint.h"
#include <iostream>

IVPStatus ExplicitMidpoint::integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept{
    Eigen::VectorXd k1 = f(t, ic);
    Eigen::VectorXd k2 = f(t+0.5*dt, ic+0.5*dt*k1);
    ic += dt*k2;

    if (Eigen::isfinite(ic.array()).all()) return IVPStatus::Success;
    return IVPStatus::FailureToSolve;
}