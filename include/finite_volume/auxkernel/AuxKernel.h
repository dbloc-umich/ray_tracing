// An auxiliary kernel computes a quantity of interest given a state vector within a cell. That quantity can be a scalar or vector.
#ifndef AUXKERNEL_H
#define AUXKERNEL_H

#include "StateMesh.h"
template<int N = 1>
class AuxKernel{
    public:
    static_assert(N >= 1, "N must be a non-zero number.");
    using RangeType = std::conditional_t<N == 1, double, Eigen::Vector<double, N>>;
    
    AuxKernel() {}
    virtual ~AuxKernel() = default;
    // Best if the map is returned by StateMesh::stateMap
    virtual RangeType computeValue(const std::map<std::string, double>& u) const = 0;
};

#endif