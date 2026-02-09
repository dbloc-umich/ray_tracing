// Structured mesh storage class for finite-volume problems. It only stores the locations of the cell boundaries on each axis.
#ifndef FV_SPATIAL_MESH_H
#define FV_SPATIAL_MESH_H

#include "UnitVector.h"

class FVSpatialMesh{
    public:
    FVSpatialMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z);
    FVSpatialMesh(const std::array<Eigen::ArrayXd, 3>& axes);
    virtual ~FVSpatialMesh() = default;

    Eigen::Index size(Eigen::Index varID) const noexcept{ return _axes[varID].size(); }
    const Eigen::ArrayXd& axis(Eigen::Index varID) const noexcept{ return _axes[varID]; }

    virtual bool isOrthogonal() const noexcept = 0;
    virtual Eigen::Vector3d cartesian(double, double, double) const noexcept = 0; // converts current coordinates to cartesian
    virtual Eigen::Vector3d coordinates(double, double, double) const noexcept = 0; // converts cartesian to current coordinates
    virtual Eigen::Vector3d basisVector(Eigen::Index varID, double, double, double) const noexcept = 0;
    // normal vector should always point towards the direction where the variable is increasing
    virtual UnitVector3d normalVector(Eigen::Index varID, double, double, double) const noexcept = 0;

    /* surfID is an integer from 0 to 5 corresponding to the following surfaces: r0, r1, mu0, mu1, ps0, ps1 */
    virtual double area(Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const noexcept = 0;
    virtual double volume(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept = 0;

    protected:
    std::array<Eigen::ArrayXd, 3> _axes; // stores the boundary cell values on each axis
};

#endif