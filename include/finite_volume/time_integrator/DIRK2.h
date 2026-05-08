#ifndef DIRK2_H
#define DIRK2_H

#include "ImplicitTimeIntegrator.h"

class DIRK2: public ImplicitTimeIntegrator{
    public:
    using ImplicitTimeIntegrator::ImplicitTimeIntegrator;
    IVPStatus integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const override;

    protected:
    inline static const double _x = 1.0 - 1.0/std::sqrt(2);
};

#endif