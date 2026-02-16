#ifndef EXPLICIT_MIDPOINT_H
#define EXPLICIT_MIDPOINT_H

#include "FVTimeIntegrator.h"

class FVExplicitMidpoint: public FVTimeIntegrator{
    public:
    IVPStatus integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept override;
};

#endif