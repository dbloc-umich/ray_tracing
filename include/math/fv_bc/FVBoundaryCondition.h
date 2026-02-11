// This class represents a boundary condition 
#ifndef FV_BOUNDARY_CONDITION_H
#define FV_BOUNDARY_CONDITION_H

#include "Material.h"
#include "Eigen/Dense"

class FVBoundaryCondition{
    public:
    FVBoundaryCondition(std::shared_ptr<Material> mat, PropVariable var, Eigen::Index surfID): _mat(mat), _var(var), _surfID(surfID) {}
    virtual ~FVBoundaryCondition() = default;
    // Computes the gradient dotted with the outward normal on the given surface
    virtual double computeFlux(double u, const Eigen::Vector3d& r) const = 0;

    protected:
    std::shared_ptr<Material> _mat;
    PropVariable _var;
    Eigen::Index _surfID;
};

#endif