/**
 * This is the full 3D mesh that stores the states at the center of each cell
 * The underlying state vector is an M-by-N matrix, where M = number of physical quantities and N = number of cells
 * Eigen stores the variable in column-major order by default, so calling cell(i) to get all the states at a single location is cache-advantaged
**/
#ifndef FV_STATE_MESH_H
#define FV_STATE_MESH_H

#include "Eigen/Dense"
#include "Material.h"

class FVMultivariateBC;
class FVSpatialMesh;
class FVVariable;

using FVCell = Eigen::MatrixXd::ColXpr;
using FVConstCell = Eigen::MatrixXd::ConstColXpr;
using FVState = Eigen::MatrixXd::RowXpr;
using FVConstState = Eigen::MatrixXd::ConstRowXpr;

class FVStateMesh{
    public:
    // Constructs with the entire state matrix
    FVStateMesh(std::shared_ptr<FVSpatialMesh>, const Eigen::MatrixXd&,
                const std::vector<std::shared_ptr<FVMultivariateBC>>&);
    FVStateMesh(std::shared_ptr<FVSpatialMesh>, Eigen::MatrixXd&&,
                const std::vector<std::shared_ptr<FVMultivariateBC>>&);
    // Constructs with a uniform initial state 
    FVStateMesh(std::shared_ptr<FVSpatialMesh>, const Eigen::VectorXd& u0,
                const std::vector<std::shared_ptr<FVMultivariateBC>>&);
    
    Eigen::Index size() const noexcept{ return _stateMesh.size(); }
    Eigen::Index cellCount() const noexcept{ return _stateMesh.cols(); }
    Eigen::Index stateCount() const noexcept{ return _stateMesh.rows(); }

    // Accessors and modifiers
    // Get variable name or index
    std::string& stateName(Eigen::Index s) noexcept{ return _stateName[s]; }
    const std::string& stateName(Eigen::Index s) const noexcept{ return _stateName[s]; }
    std::size_t stateID(const std::string& stateName) const noexcept;

    // Reference to a single value
    double& operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    const double& operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    
    // View to a cell
    FVCell operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    FVConstCell operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    FVCell cell(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{ return this->operator()(i,j,k); }
    FVConstCell cell(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{ return this->operator()(i,j,k); }

    // View to a state
    FVState operator()(Eigen::Index s) noexcept{ return _stateMesh.row(s); }
    FVConstState operator()(Eigen::Index s) const noexcept{ return _stateMesh.row(s); }
    FVState state(Eigen::Index s) noexcept{ return this->operator()(s); }
    FVConstState state(Eigen::Index s) const noexcept{ return this->operator()(s); }
    
    // Reference/View to the entire matrix
    Eigen::MatrixXd& matrix() noexcept{ return _stateMesh; }
    const Eigen::MatrixXd& matrix() const noexcept{ return _stateMesh; }
    auto flattened() noexcept{ return Eigen::Map<Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto flattened() const noexcept{ return Eigen::Map<const Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto array() noexcept{ return _stateMesh.array(); }
    auto array() const noexcept{ return _stateMesh.array(); }
    
    // Pointer to the underlying spatial mesh
    const std::shared_ptr<FVSpatialMesh> mesh() const noexcept{ return _spatialMesh; }

    // Boundary conditions
    std::shared_ptr<FVMultivariateBC> bc(Eigen::Index surfID) const noexcept{ return _bc[surfID]; }
    void setBC(Eigen::Index surfID, std::shared_ptr<FVMultivariateBC> newBC) noexcept{ _bc[surfID] = newBC; }

    // Get a map of material property variables
    Material::PropVars matProp(std::shared_ptr<Material> mat, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;

    protected:
    std::shared_ptr<FVSpatialMesh> _spatialMesh;
    std::vector<std::string> _stateName;
    Eigen::MatrixXd _stateMesh;
    std::vector<std::shared_ptr<FVMultivariateBC>> _bc;
};

#endif