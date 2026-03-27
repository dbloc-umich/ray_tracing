// Defines a Robin BC - i.e. a(r)*u + b(r)*grad(u) dot n = f(r) on the boundary
#ifndef ROBIN_BC_H
#define ROBIN_BC_H

#include "BoundaryCondition.h"

class Material;
class SpatialMesh;
class RobinBC: public BoundaryCondition{
    public:
    RobinBC(const std::function<double(const Eigen::Vector3d&)>& f,
              const std::function<double(const Eigen::Vector3d&)>& a,
              const std::function<double(const Eigen::Vector3d&)>& b,
              std::shared_ptr<SpatialMesh> mesh, Eigen::Index surfID,
              std::shared_ptr<Material> mat, std::string var, std::string prop="none");
    double computeFlux(double u, const Eigen::Vector3d& r) const override;
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _f, _a, _b;
    std::shared_ptr<SpatialMesh> _mesh;
    std::shared_ptr<Material> _mat;
    std::string _var;
    std::string _prop; // used to scale the flux if needed
};

#endif