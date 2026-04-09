#ifndef COLLISIONAL_IONIZAATION_KERNEL
#define COLLISIONAL_IONIZAATION_KERNEL

#include "Kernel.h"
class CollisionalIonizationKernel: public Kernel{
    public:
    CollisionalIonizationKernel(std::shared_ptr<Material> mat): Kernel(mat) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;
};

#endif