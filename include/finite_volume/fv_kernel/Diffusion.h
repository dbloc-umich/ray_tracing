// Implements the viscous flux term in finite volume

#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "Kernel.h"

class Diffusion: public Kernel{
    public:
    Diffusion(std::shared_ptr<Material> mat, PropVariable var, const Eigen::VectorXi& s, Prop prop=Prop::none);
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;

    protected:
    Prop _prop; // diffusion coefficient
};
#endif