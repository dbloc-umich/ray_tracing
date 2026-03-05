// Defines a Neumann BC - k(u) grad(u) dot n = f(r) on the boundary
#ifndef NEUMANN_BC_H
#define NEUMANN_BC_H

#include "UnivariateBC.h"

class SpatialMesh;
class NeumannBC: public UnivariateBC{
    public:
    NeumannBC(const std::function<double(const Eigen::Vector3d&)>& f,
                  std::shared_ptr<SpatialMesh> mesh, Eigen::Index surfID,
                  std::shared_ptr<Material> mat, PropVariable var, Prop prop=Prop::none);
    double computeFlux(double u, const Eigen::Vector3d& r) const override;
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _f;
    std::shared_ptr<SpatialMesh> _mesh;
    Prop _prop; // used to scale the flux if needed
};

#endif