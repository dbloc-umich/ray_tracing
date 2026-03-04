// Defines a Neumann BC - k(u) grad(u) dot n = f(r) on the boundary
#ifndef FV_NEUMANN_BC_H
#define FV_NEUMANN_BC_H

#include "FVUnivariateBC.h"

class FVSpatialMesh;
class FVNeumannBC: public FVUnivariateBC{
    public:
    FVNeumannBC(const std::function<double(const Eigen::Vector3d&)>& f,
                  std::shared_ptr<FVSpatialMesh> mesh, Eigen::Index surfID,
                  std::shared_ptr<Material> mat, PropVariable var, Prop prop=Prop::none);
    double computeFlux(double u, const Eigen::Vector3d& r) const override;
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _f;
    std::shared_ptr<FVSpatialMesh> _mesh;
    Prop _prop; // used to scale the flux if needed
};

#endif