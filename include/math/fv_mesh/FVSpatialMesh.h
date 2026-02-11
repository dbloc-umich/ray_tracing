// Structured mesh storage class for finite-volume problems. It only stores the locations of the cell boundaries on each axis.
#ifndef FV_SPATIAL_MESH_H
#define FV_SPATIAL_MESH_H

#include "UnitVector.h"

class FVSpatialMesh{
    public:
    FVSpatialMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z);
    FVSpatialMesh(const std::array<Eigen::ArrayXd, 3>& axes);
    virtual ~FVSpatialMesh() = default;

    const Eigen::ArrayXd& axis(Eigen::Index varID) const noexcept{ return _axes[varID]; }
    Eigen::Index axisSize(Eigen::Index varID) const noexcept{ return _axes[varID].size(); }

    virtual bool isOrthogonal() const noexcept = 0;
    virtual Eigen::Vector3d cartesian(double, double, double) const noexcept = 0; // converts current coordinates to cartesian
    Eigen::Vector3d cartesian(const Eigen::Vector3d& x) const noexcept{ return cartesian(x[0], x[1], x[2]); }
    virtual Eigen::Vector3d coordinates(double, double, double) const noexcept = 0; // converts cartesian to current coordinates
    Eigen::Vector3d coordinates(const Eigen::Vector3d& x) const noexcept{ return coordinates(x[0], x[1], x[2]); }
    virtual Eigen::Vector3d basisVector(Eigen::Index varID, double, double, double) const noexcept = 0;
    // normal vector should always point towards the direction where the variable is increasing
    virtual UnitVector3d normalVector(Eigen::Index varID, double, double, double) const noexcept = 0;
    // Scale factor/Lamé coefficient multiplied in front of the partial derivative in order to obtain the gradient
    // Use this only with orthogonal coordinate systems
    virtual double gradientFactor(Eigen::Index varID, double, double, double) const noexcept = 0;

    // surfID is an integer from 0 to 5 corresponding to the following surfaces: r0, r1, mu0, mu1, ps0, ps1
    // Grad dot normal vector, dudq is the central difference approximation of the derivative across surface #surfID at cell i,j,k
    virtual double gradientDotN(double dudq, Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const = 0;
    virtual double area(Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const noexcept = 0;
    virtual double volume(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept = 0;
    // The boundary of a ghost cell extending past cell i,j,k such that cells i,j,k and i+1,j,k are equivolume 
    virtual double ghostLength(Eigen::Index surfID) const noexcept = 0;

    protected:
    std::array<Eigen::ArrayXd, 3> _axes; // stores the boundary cell values on each axis
};

#endif