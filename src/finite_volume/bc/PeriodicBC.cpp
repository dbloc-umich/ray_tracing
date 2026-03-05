#include "PeriodicBC.h"

Eigen::VectorXd PeriodicBC::computeFlux(const ConstCell& u, const Eigen::Vector3d& r) const{
    return Eigen::VectorXd::Zero(u.rows());
}

Eigen::VectorXd PeriodicBC::computeBoundaryStates(const ConstCell& u, const Eigen::Vector3d& r) const{
    // These are dummy zeros and not physical values
    return Eigen::VectorXd::Zero(u.rows());
}