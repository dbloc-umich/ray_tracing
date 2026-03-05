#include "DegenerateBC.h"

Eigen::VectorXd DegenerateBC::computeFlux(const ConstCell& u, const Eigen::Vector3d& r) const{
    return Eigen::VectorXd::Zero(u.rows());
}

Eigen::VectorXd DegenerateBC::computeBoundaryStates(const ConstCell& u, const Eigen::Vector3d& r) const{
    // These are dummy zeros and not physical values
    return Eigen::VectorXd::Zero(u.rows());
}