// This class represents a single term in the residual (right-hand side) of the finite volume formulation.
// Derived classes to be implemented are Advection, DiffusionKernel, and LaserSourceKernel.

#ifndef RESIDUAL_H
#define RESIDUAL_H

#include "Eigen/Dense"

class Material;
class StateMesh;
class Kernel{
    public:
    Kernel(std::shared_ptr<Material> mat): _mat(mat) { assert(_mat); }
    virtual ~Kernel() = default;
    virtual Eigen::MatrixXd computeResidual(const StateMesh& u) const = 0;
    virtual Eigen::VectorXi stateID(const StateMesh& u) const noexcept = 0; // indices of affected variables

    protected:
    std::shared_ptr<Material> _mat;
};

#endif