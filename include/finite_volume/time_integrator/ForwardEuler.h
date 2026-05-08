#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H

#include "ExplicitTimeIntegrator.h"

class ForwardEuler: public ExplicitTimeIntegrator{
    public:
    IVPStatus integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const override;
};

#endif