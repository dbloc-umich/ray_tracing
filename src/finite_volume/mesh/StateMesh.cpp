#include "StateMesh.h"
#include "SpatialMesh.h"
#include <algorithm>
#include <unordered_set>

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh):
    _spatialMesh(spatialMesh)
{}

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, const Eigen::VectorXd& u0):
    _spatialMesh(spatialMesh)
{
    Eigen::Index M = _spatialMesh->axisSize(0)-1;
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    _stateMesh.resize(u0.size(), M*N*K);
    for (Eigen::Index i = 0; i < _stateMesh.cols(); i++) _stateMesh.col(i) = u0;

    _bc.resize(u0.size(), Eigen::NoChange);
    _bc.setConstant(nullptr);
}

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, const Eigen::MatrixXd& array):
    _spatialMesh(spatialMesh),
    _stateMesh(array)
{
    _bc.resize(array.rows(), Eigen::NoChange);
    _bc.setConstant(nullptr);
}

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, Eigen::MatrixXd&& array):
    _spatialMesh(spatialMesh),
    _stateMesh(std::move(array))
{
    _bc.resize(array.rows(), Eigen::NoChange);
    _bc.setConstant(nullptr);
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

void StateMesh::setStateName(const std::vector<std::string>& names){
    if (int(names.size()) != stateCount()) throw std::invalid_argument("ERROR: Number of names must be equal to number of states.");
    std::unordered_set<std::string> seen(names.cbegin(), names.cend());
    if (seen.size() < names.size()) throw std::invalid_argument("ERROR: Duplicate state names detected.");
    _stateName = names;
}

void StateMesh::setStateName(Eigen::Index s, std::string name){
    if (std::find(_stateName.cbegin(), _stateName.cend(), name) == _stateName.cend()) _stateName[s] = std::move(name);
    else throw std::invalid_argument("ERROR: Duplicate state names detected.");
}

void StateMesh::addVariable(std::string name, double u0){
    if (std::find(_stateName.cbegin(), _stateName.cend(), name) == _stateName.cend()){
        Eigen::Index Ns = stateCount();
        _stateName.push_back(std::move(name));
        if (_stateMesh.size() == 0){
            Eigen::Index M = _spatialMesh->axisSize(0)-1;
            Eigen::Index N = _spatialMesh->axisSize(1)-1;
            Eigen::Index K = _spatialMesh->axisSize(2)-1;
            _stateMesh.resize(1, M*N*K);
        } else _stateMesh.conservativeResize(Ns+1, Eigen::NoChange);
        _stateMesh.row(Ns).fill(u0);

        _bc.conservativeResize(Ns+1, Eigen::NoChange);
        _bc.row(Ns).fill(nullptr);
    } else throw std::invalid_argument("ERROR: Duplicate state names detected.");
}

std::map<std::string, double> StateMesh::matProp(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    std::map<std::string, double> props;
    for (Eigen::Index s = 0; s < stateCount(); s++) props[_stateName[s]] = operator()(s,i,j,k);
    return props;
}

std::map<std::string, double> StateMesh::matProp(const ConstCell& cell) const noexcept{
    std::map<std::string, double> props;
    for (Eigen::Index s = 0; s < stateCount(); s++) props[_stateName[s]] = cell[s];
    return props;
}