#include "DirichletBC.h"
#include "Material.h"
#include "SpatialMesh.h"

DirichletBC::DirichletBC(const std::function<double(const Eigen::Vector3d&)>& f,
                         std::shared_ptr<SpatialMesh> mesh, Eigen::Index surfID,
                         std::shared_ptr<Material> mat, std::string var, std::string prop):
    BoundaryCondition(surfID),
    _ub(f),
    _mesh(mesh),
    _mat(mat),
    _var(std::move(var)),
    _prop(std::move(prop))
{}

double DirichletBC::computeFlux(double u, const Eigen::Vector3d& r) const{
    // r has to be in the same coordinate system as the problem
    double ub = _ub(r);
    double ug = 2*ub - u; // ghost cell value;
    double scale, scaleg;
    
    if (_prop == "none"){
        scale = 1.0;
        scaleg = 1.0;
    } else{
        std::map<std::string, double> vars{{_var, u}};
        scale = _mat->computeProperty(_prop, vars);
        vars.at(_var) = ug;
        scaleg = _mat->computeProperty(_prop, vars);
    }

    double du = scaleg*ug - scale*u;
    double dr = (_mesh->ghostLength(_surfID) - _mesh->axis(_surfID/2)[_mesh->axis(_surfID/2).size()-2])/2;
    Eigen::Index i = std::upper_bound(_mesh->axis(0).cbegin(), _mesh->axis(0).cend(), r[0]) - _mesh->axis(0).cbegin() - 1;
    Eigen::Index j = std::upper_bound(_mesh->axis(1).cbegin(), _mesh->axis(1).cend(), r[1]) - _mesh->axis(1).cbegin() - 1;
    Eigen::Index k = std::upper_bound(_mesh->axis(2).cbegin(), _mesh->axis(2).cend(), r[2]) - _mesh->axis(2).cbegin() - 1;
    return _mesh->gradientDotN(du/dr, i, j, k, _surfID);
}

double DirichletBC::computeBoundaryState(double u, const Eigen::Vector3d& r) const{
    return _ub(r);
}
