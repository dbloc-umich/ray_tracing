#include "StateMesh.h"
#include "SpatialMesh.h"
#include <algorithm>

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, const Eigen::MatrixXd& array,
                         const std::vector<std::shared_ptr<MultivariateBC>>& bc):
    _spatialMesh(spatialMesh),
    _stateMesh(array),
    _bc(bc)
{
    assert(_bc.size() == 6);
}

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, Eigen::MatrixXd&& array,
                         const std::vector<std::shared_ptr<MultivariateBC>>& bc):
    _spatialMesh(spatialMesh),
    _stateMesh(std::move(array)),
    _bc(bc)
{
    assert(_bc.size() == 6);
}


StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, const Eigen::VectorXd& u0,
                         const std::vector<std::shared_ptr<MultivariateBC>>& bc):
    _spatialMesh(spatialMesh),
    _bc(bc)
{
    assert(_bc.size() == 6);

    Eigen::Index M = _spatialMesh->axisSize(0)-1;
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    _stateMesh.resize(u0.size(), M*N*K);
    for (Eigen::Index i = 0; i < _stateMesh.cols(); i++) _stateMesh.col(i) = u0;
}

std::size_t StateMesh::stateID(const std::string& name) const noexcept{
    auto it = std::find(_stateName.cbegin(), _stateName.cend(), name);
    return it - _stateName.cbegin();
}

double& StateMesh::operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    return _stateMesh(s, (i*N + j)*K + k);
}

const double& StateMesh::operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    return _stateMesh(s, (i*N + j)*K + k);
}

Cell StateMesh::operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    return _stateMesh.col((i*N + j)*K + k);    
}

ConstCell StateMesh::operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    return _stateMesh.col((i*N + j)*K + k);    
}

Material::PropVars StateMesh::matProp(std::shared_ptr<Material> mat, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    Material::PropVars props;
    
    Eigen::Index s = stateID("density");
    assert(s < stateCount()); // density variable has to exist
    double rho = this->operator()(s,i,j,k);
    props[PropVariable::density] = rho;

    double momentum = 0;
    for (auto var: {"x-momentum, y-momentum, z-momentum"}){
        s = stateID(var);
        if (s < stateCount()){
            double mom = this->operator()(s,i,j,k);
            momentum += mom*mom;
        }
    }
    momentum = std::sqrt(momentum);
    props[PropVariable::velocity] = momentum/rho;

    s = stateID("energy");
    if (s < stateCount()){
        double H = this->operator()(s,i,j,k)/rho;
        props[PropVariable::temperature] = mat->T_from_H(H);
    }
    return props;
}