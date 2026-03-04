#include "FVCompoundBC.h"
#include "FVUnivariateBC.h"

Eigen::VectorXd FVCompoundBC::computeFlux(const FVConstCell& u, const Eigen::Vector3d& r) const{
    assert(u.rows() == Eigen::Index(_univarBC.size()));
    Eigen::VectorXd flux(u.rows());
    for (Eigen::Index i = 0; i < u.rows(); i++) flux[i] = _univarBC[i]->computeFlux(u[i], r);
    return flux;
}

Eigen::VectorXd FVCompoundBC::computeBoundaryStates(const FVConstCell& u, const Eigen::Vector3d& r) const{
    assert(u.rows() == Eigen::Index(_univarBC.size()));
    Eigen::VectorXd states(u.rows());
    for (Eigen::Index i = 0; i < u.rows(); i++) states[i] = _univarBC[i]->computeBoundaryState(u[i], r);
    return states;
}