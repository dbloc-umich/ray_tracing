// Defines a periodic BC - i.e. u(0) and u(Nu-1) commutes with each other. The flux is actually not computed here and must be done so in a different residual object.
#ifndef PV_PERIODIC_BC
#define PV_PERIODIC_BC

#include "FVBoundaryCondition.h"

class FVPeriodicBC: public FVBoundaryCondition{
    public:
    FVPeriodicBC(PropVariable var, Eigen::Index surfID): FVBoundaryCondition(nullptr, var, surfID) {};
    double computeFlux(double u, const Eigen::Vector3d& r) const override{ return 0.0; };
};

#endif