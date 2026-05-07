#ifndef ELECTRON_HP_ENERGY_TRANSFER_KERNEL_H
#define ELECTRON_HP_ENERGY_TRANSFER_KERNEL_H

#include "Kernel.h"
class EquationOfState;
class ElectronHPEnergyTransferKernel: public Kernel{
    public:
    ElectronHPEnergyTransferKernel(std::shared_ptr<Material> mat, std::shared_ptr<Material> eMat);
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;

    protected:
    std::shared_ptr<Material> _eMat;
};

#endif