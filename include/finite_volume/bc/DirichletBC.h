// Defines a Dirichlet BC - i.e. u = f(r) on the boundary
#ifndef DIRICHLET_BC_H
#define DIRICHLET_BC_H

#include "BoundaryCondition.h"
#include <string>

class Material;
class SpatialMesh;
class DirichletBC: public BoundaryCondition{
    public:
    DirichletBC(const std::function<double(const Eigen::Vector3d&)>& f,
                  std::shared_ptr<SpatialMesh> mesh, Eigen::Index surfID,
                  std::shared_ptr<Material> mat, std::string var, std::string prop="none");
    double computeFlux(double u, const Eigen::Vector3d& r) const override;
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override;

    protected:
    std::function<double(const Eigen::Vector3d&)> _ub;
    std::shared_ptr<SpatialMesh> _mesh;
    std::shared_ptr<Material> _mat;
    std::string _var;
    std::string _prop; // used to scale the flux if needed
};

#endif