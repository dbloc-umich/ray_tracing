// This class represents a boundary condition of a single variable
#ifndef FV_UNIVARIATE_BC_H
#define FV_UNIVARIATE_BC_H

#include "Material.h"
#include "Eigen/Dense"

class FVUnivariateBC{
    public:
    FVUnivariateBC(std::shared_ptr<Material> mat, PropVariable var, Eigen::Index surfID): _mat(mat), _var(var), _surfID(surfID) {}
    virtual ~FVUnivariateBC() = default;
    // Computes the gradient dotted with the outward normal on the given surface
    virtual double computeFlux(double u, const Eigen::Vector3d& r) const = 0;
    virtual double computeBoundaryState(double u, const Eigen::Vector3d& r) const = 0;

    protected:
    std::shared_ptr<Material> _mat;
    PropVariable _var;
    Eigen::Index _surfID;
};

#endif