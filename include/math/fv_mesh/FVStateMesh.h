// This is the full 3D mesh that stores a physical quantity at the center of each cell
#ifndef FV_STATE_MESH_H
#define FV_STATE_MESH_H

#include "Eigen/Dense"

class FVSpatialMesh;
class FVStateMesh{
    public:
    FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, const Eigen::VectorXd& array);
    FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, Eigen::VectorXd&& array);
    FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, double u0=0);
    
    Eigen::Index size(){ return _stateMesh.size(); }

    // Accessors and modifiers
    double& operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    const double& operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    Eigen::VectorXd& matrix() noexcept{ return _stateMesh; }
    const Eigen::VectorXd& matrix() const noexcept{ return _stateMesh; }
    Eigen::ArrayXd array() noexcept{ return _stateMesh.array(); }
    const Eigen::ArrayXd array() const noexcept{ return _stateMesh.array(); }
    const std::shared_ptr<FVSpatialMesh> mesh() const noexcept{ return _spatialMesh; }

    protected:
    std::shared_ptr<FVSpatialMesh> _spatialMesh;
    Eigen::VectorXd _stateMesh;
};

#endif