// Implements the viscous flux term in finite volume

#ifndef DIFFUSION_KERNEL_H
#define DIFFUSION_KERNEL_H

#include "Kernel.h"
class AuxKernel;
class DiffusionKernel: public Kernel{
    public:
    DiffusionKernel(std::shared_ptr<Material> mat, std::string var,
                    std::shared_ptr<AuxKernel>, std::string primitive="none", std::string prop="none");
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept;

    protected:
    std::string _var; // affected variable in the system
    std::shared_ptr<AuxKernel> _converter; // auxkernel to convert to a primitive variable
    std::string _primitive; // primitive variable name
    std::string _prop; // diffusion coefficient
};
#endif