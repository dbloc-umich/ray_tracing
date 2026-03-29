#ifndef ELECTRON_HP_ENERGY_TRANSFER_KERNEL_H
#define ELECTRON_HP_ENERGY_TRANSFER_KERNEL_H

#include "Kernel.h"
class EquationOfState;
class ElectronHPEnergyTransferKernel: public Kernel{
    public:
    ElectronHPEnergyTransferKernel(std::shared_ptr<Material> mat, std::shared_ptr<EquationOfState> eos);
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept;

    protected:
    std::shared_ptr<EquationOfState> _eos;
};

#endif