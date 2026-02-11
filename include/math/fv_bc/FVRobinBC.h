// Defines a Robin BC - i.e. a(r)*u + b(r)*grad(u) dot n = f(r) on the boundary
#ifndef FV_ROBIN_BC_H
#define FV_ROBIN_BC_H

#include "FVBoundaryCondition.h"

class FVSpatialMesh;
class FVRobinBC: public FVBoundaryCondition{
    public:
    FVRobinBC(const std::function<double(const Eigen::Vector3d&)>& f,
              const std::function<double(const Eigen::Vector3d&)>& a,
              const std::function<double(const Eigen::Vector3d&)>& b,
              std::shared_ptr<FVSpatialMesh> mesh, Eigen::Index surfID,
              std::shared_ptr<Material> mat, PropVariable var, Prop prop=Prop::none);
    double computeFlux(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _f, _a, _b;
    std::shared_ptr<FVSpatialMesh> _mesh;
    Prop _prop; // used to scale the flux if needed
};

#endif