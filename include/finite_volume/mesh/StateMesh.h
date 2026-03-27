/**
 * This is the full 3D mesh that stores the states at the center of each cell
 * The underlying state vector is an M-by-N matrix, where M = number of physical quantities and N = number of cells
 * Eigen stores the variable in column-major order by default, so calling cell(i) to get all the states at a single location is cache-advantaged
**/
#ifndef STATE_MESH_H
#define STATE_MESH_H

#include "Eigen/Dense"
#include "Material.h"

class BoundaryCondition;
class SpatialMesh;
class Variable;

using Cell = Eigen::MatrixXd::ColXpr;
using ConstCell = Eigen::MatrixXd::ConstColXpr;
using State = Eigen::MatrixXd::RowXpr;
using ConstState = Eigen::MatrixXd::ConstRowXpr;

class StateMesh{
    public:
    // Constructs without any initial information
    StateMesh(std::shared_ptr<SpatialMesh>);
    // Constructs with a uniform initial state 
    StateMesh(std::shared_ptr<SpatialMesh>, const Eigen::VectorXd& u0);
    // Constructs with the entire state matrix
    StateMesh(std::shared_ptr<SpatialMesh>, const Eigen::MatrixXd&);
    StateMesh(std::shared_ptr<SpatialMesh>, Eigen::MatrixXd&&);
    
    Eigen::Index size() const noexcept{ return _stateMesh.size(); }
    Eigen::Index cellCount() const noexcept{ return _stateMesh.cols(); }
    Eigen::Index stateCount() const noexcept{ return _stateMesh.rows(); }

    // Accessors and modifiers
    // Reference to a single value
    double& operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    const double& operator()(Eigen::Index s, Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    
    // View to a cell
    Cell operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept;
    ConstCell operator()(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    Cell cell(Eigen::Index i, Eigen::Index j, Eigen::Index k) noexcept{ return this->operator()(i,j,k); }
    ConstCell cell(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{ return this->operator()(i,j,k); }

    // View to a state
    State operator()(Eigen::Index s) noexcept{ return _stateMesh.row(s); }
    ConstState operator()(Eigen::Index s) const noexcept{ return _stateMesh.row(s); }
    State state(Eigen::Index s) noexcept{ return this->operator()(s); }
    ConstState state(Eigen::Index s) const noexcept{ return this->operator()(s); }
    
    // Reference/View to the entire matrix
    Eigen::MatrixXd& matrix() noexcept{ return _stateMesh; }
    const Eigen::MatrixXd& matrix() const noexcept{ return _stateMesh; }
    auto flattened() noexcept{ return Eigen::Map<Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto flattened() const noexcept{ return Eigen::Map<const Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto array() noexcept{ return _stateMesh.array(); }
    auto array() const noexcept{ return _stateMesh.array(); }
    
    // Pointer to the underlying spatial mesh
    const std::shared_ptr<SpatialMesh> mesh() const noexcept{ return _spatialMesh; }

    // Get variable name or index
    std::string& stateName(Eigen::Index s) noexcept{ return _stateName[s]; }
    const std::string& stateName(Eigen::Index s) const noexcept{ return _stateName[s]; }
    std::size_t stateID(const std::string& stateName) const noexcept;
    const std::vector<std::string>& stateName() const noexcept{ return _stateName; }
    void setStateName(const std::vector<std::string>& name);
    void setStateName(Eigen::Index s, std::string name);
    void addVariable(std::string name, double u0);

    // Boundary conditions
    std::shared_ptr<BoundaryCondition> bc(Eigen::Index varID, Eigen::Index surfID) const noexcept{ return _bc(varID, surfID); }
    void setBC(Eigen::Index varID, Eigen::Index surfID, std::shared_ptr<BoundaryCondition> bc) noexcept{ _bc(varID, surfID) = bc; }

    // Get a map of material property variables
    std::map<std::string, double> matProp(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept;
    std::map<std::string, double> matProp(const ConstCell&) const noexcept;

    protected:
    std::shared_ptr<SpatialMesh> _spatialMesh;
    std::vector<std::string> _stateName;
    Eigen::MatrixXd _stateMesh;
    Eigen::Array<std::shared_ptr<BoundaryCondition>, Eigen::Dynamic, 6> _bc;
};

#endif