/**
 * This is the full 3D mesh that stores the states at the center of each cell
 * The underlying state vector is an M-by-N matrix, where M = number of physical quantities and N = number of cells
 * Eigen stores the variable in column-major order by default, so calling cell(i) to get all the states at a single location is cache-advantaged
**/
#ifndef FV_STATE_MESH_H
#define FV_STATE_MESH_H

#include "Eigen/Dense"

class FVSpatialMesh;
class FVStateMesh{
    public:
    // Constructs with the entire state matrix
    FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, const Eigen::MatrixXd& array);
    FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, Eigen::MatrixXd&& array);
    FVStateMesh(std::shared_ptr<FVSpatialMesh> spatialMesh, Eigen::Index states=1, double u0=0);
    
    Eigen::Index size() const noexcept{ return _stateMesh.size(); }
    Eigen::Index cellCount() const noexcept{ return _stateMesh.cols(); }
    Eigen::Index stateCount() const noexcept{ return _stateMesh.rows(); }

    // Accessors and modifiers
    // Reference to a single value
    double& operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    const double& operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    // View to a cell
    Eigen::MatrixXd::ColXpr operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    Eigen::MatrixXd::ConstColXpr operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    Eigen::MatrixXd::ColXpr cell(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{ return this->operator()(i,j,k); }
    Eigen::MatrixXd::ConstColXpr cell(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{ return this->operator()(i,j,k); }
    // View to a state
    Eigen::MatrixXd::RowXpr operator()(Eigen::Index s) noexcept{ return _stateMesh.row(s); }
    Eigen::MatrixXd::ConstRowXpr operator()(Eigen::Index s) const noexcept{ return _stateMesh.row(s); }
    Eigen::MatrixXd::RowXpr state(Eigen::Index s) noexcept{ return this->operator()(s); }
    Eigen::MatrixXd::ConstRowXpr state(Eigen::Index s) const noexcept{ return this->operator()(s); }
    // Reference/View to the entire matrix
    Eigen::MatrixXd& matrix() noexcept{ return _stateMesh; }
    const Eigen::MatrixXd& matrix() const noexcept{ return _stateMesh; }
    auto flattened() noexcept{ return Eigen::Map<Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto flattened() const noexcept{ return Eigen::Map<const Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto array() noexcept{ return _stateMesh.array(); }
    auto array() const noexcept{ return _stateMesh.array(); }
    // Reference to the underlying spatial mesh
    const std::shared_ptr<FVSpatialMesh> mesh() const noexcept{ return _spatialMesh; }

    protected:
    std::shared_ptr<FVSpatialMesh> _spatialMesh;
    Eigen::MatrixXd _stateMesh;
};

#endif