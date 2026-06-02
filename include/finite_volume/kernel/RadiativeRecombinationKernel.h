#ifndef RADIATIVE_RECOMBINATION_KERNEL
#define RADIATIVE_RECOMBINATION_KERNEL

#include "Kernel.h"
class RadiativeRecombinationKernel: public Kernel{
    public:
    RadiativeRecombinationKernel(std::shared_ptr<Material> mat): Kernel(mat) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;
};

#endif