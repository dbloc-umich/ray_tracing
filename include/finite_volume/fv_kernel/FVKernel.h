// This class represents a single term in the residual (right-hand side) of the finite volume formulation.
// Derived classes to be implemented are FVAdvection, FVDiffusion, and FVLaserSource.

#ifndef FV_RESIDUAL_H
#define FV_RESIDUAL_H

#include "Material.h"
#include "Eigen/Dense"

class FVStateMesh;
class FVKernel{
    public:
    FVKernel(std::shared_ptr<Material> mat, PropVariable var, Eigen::Index s): _mat(mat), _var(var), _s(s) {}
    virtual ~FVKernel() = default;
    virtual Eigen::VectorXd computeResidual(const FVStateMesh& u) const = 0;

    Eigen::Index stateID() const noexcept{ return _s; }

    protected:
    std::shared_ptr<Material> _mat;
    PropVariable _var; // variable name used for material lookup
    Eigen::Index _s; // variable ID in the solver
};

#endif