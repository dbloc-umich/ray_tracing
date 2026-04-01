// Defines a degenerate BC, when the surface collapses into a point or a line
#ifndef DEGENERATE_BC_H
#define DEGENERATE_BC_H

#include "BoundaryCondition.h"

class DegenerateBC: public BoundaryCondition{
    public:
    DegenerateBC(Eigen::Index surfID): BoundaryCondition(surfID) {}
    double computeFlux(double u, const Eigen::Vector3d& r) const override{ return 0.0; }
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override{ return u; }
};

#endif