// Defines a Robin BC - i.e. a(r)*u + b(r)*grad(u) dot n = f(r) on the boundary
#ifndef ROBIN_BC_H
#define ROBIN_BC_H

#include "UnivariateBC.h"

class SpatialMesh;
class RobinBC: public UnivariateBC{
    public:
    RobinBC(const std::function<double(const Eigen::Vector3d&)>& f,
              const std::function<double(const Eigen::Vector3d&)>& a,
              const std::function<double(const Eigen::Vector3d&)>& b,
              std::shared_ptr<SpatialMesh> mesh, Eigen::Index surfID,
              std::shared_ptr<Material> mat, PropVariable var, Prop prop=Prop::none);
    double computeFlux(double u, const Eigen::Vector3d& r) const override;
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _f, _a, _b;
    std::shared_ptr<SpatialMesh> _mesh;
    Prop _prop; // used to scale the flux if needed
};

#endif