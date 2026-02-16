// Implements the viscous flux term in finite volume

#ifndef FV_DIFFUSION_H
#define FV_DIFFUSION_H

#include "FVResidual.h"

class FVBoundaryCondition;
class FVDiffusion: public FVResidual{
    public:
    FVDiffusion(const std::vector<std::shared_ptr<FVBoundaryCondition>>&,
                std::shared_ptr<Material> mat, PropVariable var, Prop prop=Prop::none);
    Eigen::VectorXd computeResidual(const FVStateMesh& u) const override;

    protected:
    std::vector<std::shared_ptr<FVBoundaryCondition>> _bc;
    Prop _prop; // diffusion coefficient
    std::vector<bool> _periodic;
};
#endif