#ifndef THREE_BODY_RECOMBINATION_KERNEL
#define THREE_BODY_RECOMBINATION_KERNEL

#include "Kernel.h"
class ThreeBodyRecombinationKernel: public Kernel{
    public:
    ThreeBodyRecombinationKernel(std::shared_ptr<Material> mat): Kernel(mat) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;
};

#endif