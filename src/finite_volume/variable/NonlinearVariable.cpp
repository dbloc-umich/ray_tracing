// #include "NonlinearVariable.h"
// #include "Kernel.h"
// #include "StateMesh.h"

// void NonlinearVariable::addKernel(std::unique_ptr<Kernel> kernel){
//     if (kernel) _kernels.push_back(std::move(kernel));
// }

// Eigen::VectorXd NonlinearVariable::computeResidual(const StateMesh& u) const{
//     Eigen::Index N = u.cellCount();
//     Eigen::VectorXd R = Eigen::VectorXd::Zero(N);
//     for (const auto& kernel: _kernels) R += kernel->computeResidual(u);
//     return R;
// }