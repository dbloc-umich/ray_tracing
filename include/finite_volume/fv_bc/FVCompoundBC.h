// Mutivariate boundary condition built from multiple univariable boundary conditions
#ifndef FV_COMPOUND_BC_H
#define FV_COMPOUND_BC_H

#include "FVMultivariateBC.h"

class FVUnivariateBC;
class FVCompoundBC: public FVMultivariateBC{
    public:
    FVCompoundBC(const std::vector<std::shared_ptr<FVUnivariateBC>>& univarBC): _univarBC(univarBC) {}
    Eigen::VectorXd computeFlux(const FVConstCell& u, const Eigen::Vector3d& r) const override;
    Eigen::VectorXd computeBoundaryStates(const FVConstCell& u, const Eigen::Vector3d& r) const override;

    protected:
    std::vector<std::shared_ptr<FVUnivariateBC>> _univarBC;
};

#endif