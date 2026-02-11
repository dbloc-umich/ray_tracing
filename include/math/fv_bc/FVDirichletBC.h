// Defines a Dirichlet BC - i.e. u = f(r) on the boundary
#ifndef FV_DIRICHLET_BC_H
#define FV_DIRICHLET_BC_H

#include "FVBoundaryCondition.h"

class FVSpatialMesh;
class FVDirichletBC: public FVBoundaryCondition{
    public:
    FVDirichletBC(const std::function<double(const Eigen::Vector3d&)>& f,
                  std::shared_ptr<FVSpatialMesh> mesh, Eigen::Index surfID,
                  std::shared_ptr<Material> mat, PropVariable var, Prop prop=Prop::none);
    double computeFlux(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _ub;
    std::shared_ptr<FVSpatialMesh> _mesh;
    Prop _prop; // used to scale the flux if needed
};

#endif