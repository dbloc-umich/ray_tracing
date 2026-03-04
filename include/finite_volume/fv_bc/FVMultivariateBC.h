// This class represents a boundary condition of multiple variables (i.e. the entire cell at the same time)
#ifndef FV_MULTIVARIATE_BC_H
#define FV_MULTIVARIATE_BC_H

#include "FVStateMesh.h"
class FVMultivariateBC{
    public:
    virtual ~FVMultivariateBC() = default;
    virtual Eigen::VectorXd computeFlux(const FVConstCell& u, const Eigen::Vector3d& r) const = 0;
    virtual Eigen::VectorXd computeBoundaryStates(const FVConstCell& u, const Eigen::Vector3d& r) const = 0;
    virtual bool isPeriodicBC() const noexcept{ return false; }
};

#endif