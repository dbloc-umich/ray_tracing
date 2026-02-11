#include "FVStateMesh.h"
#include "FVSpatialMesh.h"

FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, const Eigen::ArrayXd& array):
    _spatialMesh(spatialMesh),
    _stateMesh(array)
{}

FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, Eigen::ArrayXd&& array):
    _spatialMesh(spatialMesh),
    _stateMesh(std::move(array))
{}


FVStateMesh::FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, double ic):
    _spatialMesh(spatialMesh)
{
    Eigen::Index M = _spatialMesh->axisSize(0)-1;
    Eigen::Index N = _spatialMesh->axisSize(1)-1;
    Eigen::Index K = _spatialMesh->axisSize(2)-1;
    _stateMesh = Eigen::ArrayXd::Constant(M*N*K, ic);
}

// FVStateMesh::FVStateMesh(const Eigen::ArrayXd& array):
//     _spatialMesh(nullptr),
//     _stateMesh(array)
// {}

// FVStateMesh::FVStateMesh(Eigen::ArrayXd&& array):
//     _spatialMesh(nullptr),
//     _stateMesh(std::move(array))
// {}

// FVStateMesh::FVStateMesh(std::size_t sz, double ic):
//     _spatialMesh(nullptr),
//     _stateMesh(Eigen::ArrayXd::Constant(sz, ic))
// {}

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