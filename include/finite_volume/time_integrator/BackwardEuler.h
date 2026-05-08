#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H

#include "ImplicitTimeIntegrator.h"

class BackwardEuler: public ImplicitTimeIntegrator{
    public:
    using ImplicitTimeIntegrator::ImplicitTimeIntegrator;
    IVPStatus integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const override;
};

#endif