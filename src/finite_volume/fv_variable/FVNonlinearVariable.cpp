#include "FVNonlinearVariable.h"
#include "FVKernel.h"
#include "FVStateMesh.h"

void FVNonlinearVariable::addKernel(std::unique_ptr<FVKernel> kernel){
    if (kernel) _kernels.push_back(std::move(kernel));
}

Eigen::VectorXd FVNonlinearVariable::computeResidual(const FVStateMesh& u) const{
    Eigen::Index N = u.cellCount();
    Eigen::VectorXd R = Eigen::VectorXd::Zero(N);
    for (const auto& kernel: _kernels) R += kernel->computeResidual(u);
    return R;
}