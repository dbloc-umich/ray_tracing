// Defines a degenerate BC, when the surface collapses into a point or a line
#ifndef DEGENERATE_BC_H
#define DEGENERATE_BC_H

#include "MultivariateBC.h"

class DegenerateBC: public MultivariateBC{
    public:
    virtual Eigen::VectorXd computeFlux(const ConstCell& u, const Eigen::Vector3d& r) const override;
    virtual Eigen::VectorXd computeBoundaryStates(const ConstCell& u, const Eigen::Vector3d& r) const override;
};

#endif