// Implements the viscous flux term in finite volume

#ifndef FV_DIFFUSION_H
#define FV_DIFFUSION_H

#include "FVResidual.h"

class FVDiffusion: public FVResidual{
    public:
    FVDiffusion(std::shared_ptr<Material> mat, PropVariable var): FVResidual(mat, var) {}
    FVStateMesh computeResidual(const FVStateMesh& u) const override;
};
#endif