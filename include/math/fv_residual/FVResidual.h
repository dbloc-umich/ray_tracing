// This class represents a single term in the residual (right-hand side) of the finite volume formulation.
// Derived classes to be implemented are FVAdvection, FVDiffusion, and FVLaserSource.

#ifndef FV_RESIDUAL_H
#define FV_RESIDUAL_H

#include "Material.h"
#include "Eigen/Dense"

class FVStateMesh;
class FVResidual{
    public:
    FVResidual(std::shared_ptr<Material> mat, PropVariable var): _mat(mat), _var(var) {}
    virtual ~FVResidual() = default;
    virtual Eigen::VectorXd computeResidual(const FVStateMesh& u) const = 0;

    protected:
    std::shared_ptr<Material> _mat;
    PropVariable _var;
};

#endif