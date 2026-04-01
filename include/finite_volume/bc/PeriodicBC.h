// Defines a periodic BC - i.e. u(0) and u(Nu-1) commutes with each other. The flux is actually not computed here and must be done so in an Kernel object.
#ifndef PERIODIC_BC
#define PERIODIC_BC

#include "BoundaryCondition.h"

class PeriodicBC: public BoundaryCondition{
    public:
    PeriodicBC(Eigen::Index surfID): BoundaryCondition(surfID) {}
    double computeFlux(double u, const Eigen::Vector3d& r) const override{ return 0.0; }
    double computeBoundaryState(double u, const Eigen::Vector3d& r) const override{ return u; }
    bool isPeriodic() const noexcept override{ return true; }
};

#endif