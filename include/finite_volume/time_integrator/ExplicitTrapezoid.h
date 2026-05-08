#ifndef EXPLICIT_TRAPEZOID_H
#define EXPLICIT_TRAPEZOID_H

#include "ExplicitTimeIntegrator.h"

class ExplicitTrapezoid: public ExplicitTimeIntegrator{
    public:
    IVPStatus integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const override;
};

#endif