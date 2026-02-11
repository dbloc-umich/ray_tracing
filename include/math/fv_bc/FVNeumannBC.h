// Defines a Neumann BC - k(u) grad(u) dot n = f(r) on the boundary
#ifndef FV_NEUMANN_BC_H
#define FV_NEUMANN_BC_H

#include "FVBoundaryCondition.h"

class FVNeumannBC: public FVBoundaryCondition{
    public:
    FVNeumannBC(const std::function<double(const Eigen::Vector3d&)>& f, PropVariable var, Eigen::Index surfID);
    double computeFlux(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _f;
};

#endif