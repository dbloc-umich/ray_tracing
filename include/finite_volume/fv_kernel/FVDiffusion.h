// Implements the viscous flux term in finite volume

#ifndef FV_DIFFUSION_H
#define FV_DIFFUSION_H

#include "FVKernel.h"

class FVDiffusion: public FVKernel{
    public:
    FVDiffusion(std::shared_ptr<Material> mat, PropVariable var, Eigen::Index s, Prop prop=Prop::none);
    Eigen::VectorXd computeResidual(const FVStateMesh& u) const override;

    protected:
    Prop _prop; // diffusion coefficient
};
#endif