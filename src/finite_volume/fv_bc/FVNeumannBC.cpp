#include "FVNeumannBC.h"
#include "FVSpatialMesh.h"

FVNeumannBC::FVNeumannBC(const std::function<double(const Eigen::Vector3d&)>& f,
                         std::shared_ptr<FVSpatialMesh> mesh, Eigen::Index surfID,
                         std::shared_ptr<Material> mat, PropVariable var, Prop prop):
    FVUnivariateBC(mat, var, surfID),
    _f(f),
    _mesh(mesh),
    _prop(prop)
{}

double FVNeumannBC::computeFlux(double u, const Eigen::Vector3d& r) const{
    // r has to be in the same coordinate system as the problem
    return _f(r);
}

double FVNeumannBC::computeBoundaryState(double u, const Eigen::Vector3d& r) const{
    return 0.0;
    // double flux = _f(r);
    // Eigen::Index i = std::upper_bound(_mesh->axis(0).cbegin(), _mesh->axis(0).cend(), r[0]) - _mesh->axis(0).cbegin() - 1;
    // Eigen::Index j = std::upper_bound(_mesh->axis(1).cbegin(), _mesh->axis(1).cend(), r[1]) - _mesh->axis(1).cbegin() - 1;
    // Eigen::Index k = std::upper_bound(_mesh->axis(2).cbegin(), _mesh->axis(2).cend(), r[2]) - _mesh->axis(2).cbegin() - 1;
    
    // double dr = (_mesh->ghostLength(_surfID) - _mesh->axis(_surfID/2)[_mesh->axis(_surfID/2).size()-2])/2;
    // if (_surfID )
    // switch (_prop){
    //     case Prop::none: {
    //         scale = 1.0;
    //         scaleg = 1.0;
    //         break;
    //     }
    //     default: {
    //         Material::PropVars vars{{_var, u}};
    //         scale = _mat->computeProperty(_prop, vars);
    //         vars.at(_var) = ug;
    //         scaleg = _mat->computeProperty(_prop, vars);
    //         break;
    //     }
    // }

    // double du = scaleg*ug - scale*u;

    // return _mesh->gradientDotN(du/dr, i, j, k, _surfID);
}