#include "FVNeumannBC.h"

FVNeumannBC::FVNeumannBC(const std::function<double(const Eigen::Vector3d&)>& f, PropVariable var, Eigen::Index surfID):
    FVBoundaryCondition(nullptr, var, surfID),
    _f(f)
{}

double FVNeumannBC::computeFlux(double u, const Eigen::Vector3d& r) const{
    // r has to be in the same coordinate system as the problem
    return _f(r);
}