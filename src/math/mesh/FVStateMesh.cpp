#include "FVStateMesh.h"
#include "FVSpatialMesh.h"

FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, double ic):
    _spatialMesh(spatialMesh)
{
    Eigen::Index M = _spatialMesh->size(0);
    Eigen::Index N = _spatialMesh->size(1);
    Eigen::Index K = _spatialMesh->size(2);
    _stateMesh = Eigen::ArrayXd::Constant((M-1)*(N-1)*(K-1), ic);
}

double& FVStateMesh::operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{
    Eigen::Index M = _spatialMesh->size(0);
    Eigen::Index N = _spatialMesh->size(1);
    return _stateMesh[(i*N + j)*M + k];
}

const double& FVStateMesh::operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    Eigen::Index M = _spatialMesh->size(0);
    Eigen::Index N = _spatialMesh->size(1);
    return _stateMesh[(i*N + j)*M + k];
}