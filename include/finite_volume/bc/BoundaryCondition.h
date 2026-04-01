// This class represents a boundary condition of a single variable
#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Eigen/Dense"

class BoundaryCondition{
    public:
    BoundaryCondition(Eigen::Index surfID): _surfID(surfID) {}
    virtual ~BoundaryCondition() = default;
    // Computes the gradient dotted with the outward normal on the given surface
    virtual double computeFlux(double u, const Eigen::Vector3d& r) const = 0;
    virtual double computeBoundaryState(double u, const Eigen::Vector3d& r) const = 0;
    virtual bool isPeriodic() const noexcept{ return false; }

    Eigen::Index surfID() const noexcept{ return _surfID; }

    protected:
    Eigen::Index _surfID;
};

#endif