// Mutivariate boundary condition built from multiple univariable boundary conditions
#ifndef COMPOUND_BC_H
#define COMPOUND_BC_H

#include "MultivariateBC.h"

class UnivariateBC;
class CompoundBC: public MultivariateBC{
    public:
    CompoundBC(const std::vector<std::shared_ptr<UnivariateBC>>& univarBC): _univarBC(univarBC) {}
    Eigen::VectorXd computeFlux(const ConstCell& u, const Eigen::Vector3d& r) const override;
    Eigen::VectorXd computeBoundaryStates(const ConstCell& u, const Eigen::Vector3d& r) const override;

    protected:
    std::vector<std::shared_ptr<UnivariateBC>> _univarBC;
};

#endif