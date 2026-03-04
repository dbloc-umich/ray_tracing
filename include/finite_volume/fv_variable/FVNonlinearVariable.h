#ifndef FV_NONLINEAR_VARIABLE_H
#define FV_NONLINEAR_VARIABLE_H

#include "FVVariable.h"
#include "Eigen/Dense"
#include <memory>
#include <vector>

class FVKernel;
class FVStateMesh;
class FVNonlinearVariable: public FVVariable{
    public:
    FVNonlinearVariable(const std::string& name, unsigned short dim = 1): FVVariable(name, dim) {}
    bool aux() const noexcept override{ return false; }
    void addKernel(std::unique_ptr<FVKernel>);
    Eigen::VectorXd computeResidual(const FVStateMesh& u) const;

    protected:
    std::vector<std::unique_ptr<FVKernel>> _kernels;
};

#endif