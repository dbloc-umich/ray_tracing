#include "FVPeriodicBC.h"

Eigen::VectorXd FVPeriodicBC::computeFlux(const FVConstCell& u, const Eigen::Vector3d& r) const{
    return Eigen::VectorXd::Zero(u.rows());
}

Eigen::VectorXd FVPeriodicBC::computeBoundaryStates(const FVConstCell& u, const Eigen::Vector3d& r) const{
    // These are dummy zeros and not physical values
    return Eigen::VectorXd::Zero(u.rows());
}