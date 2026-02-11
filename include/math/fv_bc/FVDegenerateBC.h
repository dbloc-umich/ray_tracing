// Defines a degenerate BC, when the surface collapses into a point or a line
#ifndef FV_DEGENERATE_BC_H
#define FV_DEGENERATE_BC_H

#include "FVBoundaryCondition.h"

class FVDegenerateBC: public FVBoundaryCondition{
    public:
    FVDegenerateBC(PropVariable var, Eigen::Index surfID): FVBoundaryCondition(nullptr, var, surfID) {};
    double computeFlux(double, const Eigen::Vector3d&) const override{ return 0.0; };
};

#endif