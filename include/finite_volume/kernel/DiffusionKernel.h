// Implements the viscous flux term in finite volume

#ifndef DIFFUSION_KERNEL_H
#define DIFFUSION_KERNEL_H

#include "Kernel.h"
template<int N>
class AuxKernel;

template<int N>
class DiffusionKernel: public Kernel{
    public:
    static_assert(N >= 1 && N <= 3, "Diffusion is supported for up to 3 variables.");
    using RangeType = std::conditional_t<N == 1, double, Eigen::Vector<double, N>>;
    
    DiffusionKernel(std::shared_ptr<Material> mat, std::string var,
                    std::shared_ptr<AuxKernel<N>>, std::string primitive="none", std::string prop="none");
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;

    protected:
    std::string _var; // affected variable in the system
    std::shared_ptr<AuxKernel<N>> _converter; // auxkernel to convert to a primitive variable
    std::string _primitive; // primitive variable name
    std::string _prop; // diffusion coefficient
};
#endif