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

    _L.resize(u0.size());
    _L.setConstant(-std::numeric_limits<double>::infinity());
    
    _U.resize(u0.size());
    _U.setConstant(std::numeric_limits<double>::infinity());
}

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, const Eigen::MatrixXd& array):
    _spatialMesh(spatialMesh),
    _stateMesh(array)
{
    _bc.resize(array.rows(), Eigen::NoChange);
    _bc.setConstant(nullptr);

    _L.resize(array.rows());
    _L.setConstant(-std::numeric_limits<double>::infinity());
    
    _U.resize(array.rows());
    _U.setConstant(std::numeric_limits<double>::infinity());
}

StateMesh::StateMesh(std::shared_ptr<SpatialMesh> spatialMesh, Eigen::MatrixXd&& array):
    _spatialMesh(spatialMesh),
    _stateMesh(std::move(array))
{
    _bc.resize(array.rows(), Eigen::NoChange);
    _bc.setConstant(nullptr);

    _L.resize(array.rows());
    _L.setConstant(-std::numeric_limits<double>::infinity());
    
    _U.resize(array.rows());
    _U.setConstant(std::numeric_limits<double>::infinity());
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

void StateMesh::addVariable(std::string name, double u0, double L, double U){
    if (L >= U) throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
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
        _L.conservativeResize(Ns+1);
        _L[Ns] = L;
        _U.conservativeResize(Ns+1);
        _U[Ns] = U;
    } else throw std::invalid_argument("ERROR: Duplicate state names detected.");
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

void StateMesh::setLowerBound(Eigen::Index s, double L){
    if (_U[s] > L) _L[s] = L;
    else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
}

void StateMesh::setLowerBound(const Eigen::VectorXd& L){
    assert(L.size() == _L.size());
    if ((_U.array() > L.array()).all()) _L = L;
    else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
}

void StateMesh::setUpperBound(Eigen::Index s, double U){
    if (_L[s] < U) _U[s] = U;
    else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
}

void StateMesh::setUpperBound(const Eigen::VectorXd& U){
    assert(U.size() == _U.size());
    if ((_L.array() < U.array()).all()) _U = U;
    else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
}

std::map<std::string, double> StateMesh::stateMap(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    std::map<std::string, double> props;
    for (Eigen::Index s = 0; s < stateCount(); s++) props[_stateName[s]] = operator()(s,i,j,k);
    return props;
}

std::map<std::string, double> StateMesh::stateMap(const ConstCell& cell) const noexcept{
    std::map<std::string, double> props;
    for (Eigen::Index s = 0; s < stateCount(); s++) props[_stateName[s]] = cell[s];
    return props;
}

void StateMesh::updateStateMap(std::map<std::string, double>& props, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    for (Eigen::Index s = 0; s < stateCount(); s++) props[_stateName[s]] = operator()(s,i,j,k);
}

void StateMesh::updateStateMap(std::map<std::string, double>& props, const ConstCell& cell) const noexcept{
    for (Eigen::Index s = 0; s < stateCount(); s++) props[_stateName[s]] = cell[s];
}