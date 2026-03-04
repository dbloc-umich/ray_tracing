// Defines a degenerate BC, when the surface collapses into a point or a line
#ifndef FV_DEGENERATE_BC_H
#define FV_DEGENERATE_BC_H

#include "FVMultivariateBC.h"

class FVDegenerateBC: public FVMultivariateBC{
    public:
    virtual Eigen::VectorXd computeFlux(const FVConstCell& u, const Eigen::Vector3d& r) const override;
    virtual Eigen::VectorXd computeBoundaryStates(const FVConstCell& u, const Eigen::Vector3d& r) const override;
};

#endif