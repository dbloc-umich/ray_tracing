// This class represents a single term in the residual (right-hand side) of the finite volume formulation.
// Derived classes to be implemented are Advection, Diffusion, and LaserSource.

#ifndef RESIDUAL_H
#define RESIDUAL_H

#include "Material.h"
#include "Eigen/Dense"

class StateMesh;
class Kernel{
    public:
    Kernel(std::shared_ptr<Material> mat, PropVariable var, const Eigen::VectorXi& s): _mat(mat), _var(var), _s(s) {}
    virtual ~Kernel() = default;
    virtual Eigen::MatrixXd computeResidual(const StateMesh& u) const = 0;

    Eigen::VectorXi stateID() const noexcept{ return _s; }

    protected:
    std::shared_ptr<Material> _mat;
    PropVariable _var; // variable name used for material lookup
    Eigen::VectorXi _s; // variable IDs in the solver
};

#endif