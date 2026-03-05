// This class represents a boundary condition of multiple variables (i.e. the entire cell at the same time)
#ifndef MULTIVARIATE_BC_H
#define MULTIVARIATE_BC_H

#include "StateMesh.h"
class MultivariateBC{
    public:
    virtual ~MultivariateBC() = default;
    virtual Eigen::VectorXd computeFlux(const ConstCell& u, const Eigen::Vector3d& r) const = 0;
    virtual Eigen::VectorXd computeBoundaryStates(const ConstCell& u, const Eigen::Vector3d& r) const = 0;
    virtual bool isPeriodicBC() const noexcept{ return false; }
};

#endif