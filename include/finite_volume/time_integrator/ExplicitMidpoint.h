#ifndef EXPLICIT_MIDPOINT_H
#define EXPLICIT_MIDPOINT_H

#include "TimeIntegrator.h"

class ExplicitMidpoint: public TimeIntegrator{
    public:
    IVPStatus integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept override;
};

#endif