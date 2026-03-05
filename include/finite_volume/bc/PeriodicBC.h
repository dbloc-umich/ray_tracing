// Defines a periodic BC - i.e. u(0) and u(Nu-1) commutes with each other. The flux is actually not computed here and must be done so in an Kernel object.
#ifndef PV_PERIODIC_BC
#define PV_PERIODIC_BC

#include "MultivariateBC.h"

class PeriodicBC: public MultivariateBC{
    public:
    Eigen::VectorXd computeFlux(const ConstCell& u, const Eigen::Vector3d& r) const override;
    Eigen::VectorXd computeBoundaryStates(const ConstCell& u, const Eigen::Vector3d& r) const override;
    bool isPeriodicBC() const noexcept override{ return true; }
};

#endif