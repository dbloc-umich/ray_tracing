#include "FVExplicitMidpoint.h"

IVPStatus FVExplicitMidpoint::integrate(const Function& f, Eigen::ArrayXd& ic, double t, double dt) const noexcept{
    Eigen::ArrayXd k1 = f(t, ic);
    Eigen::ArrayXd k2 = f(t+0.5*dt, ic+0.5*k1);
    ic += dt*k2;

    if (Eigen::isfinite(ic).all()) return IVPStatus::Success;
    return IVPStatus::FailureToSolve;
}