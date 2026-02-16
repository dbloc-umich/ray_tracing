#include "FVStateMesh.h"
#include "FVSpatialMesh.h"

FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, const Eigen::VectorXd& array):
    _spatialMesh(spatialMesh),
    _stateMesh(array)
{}

FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, Eigen::VectorXd&& array):
    _spatialMesh(spatialMesh),
    _stateMesh(std::move(array))
{}


FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, double u0):
    _spatialMesh(spatialMesh)
{
    Eigen::Index M = _spatialMesh->axisSize(0)-1;
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    _stateMesh = Eigen::VectorXd::Constant(M*N*K, u0);
}

double& FVStateMesh::operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    return _stateMesh[(i*N + j)*K + k];
}

const double& FVStateMesh::operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    return _stateMesh[(i*N + j)*K + k];
}