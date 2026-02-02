#ifndef EXPLICIT_MIDPOINT_H
#define EXPLICIT_MIDPOINT_H

#include "ODE_IVPSolver.h"

template<int M>
class ExplicitMidpoint: public ODE_IVPSolver<M>{
    public:
    using typename ODESolver<M>::RangeType;
    using typename ODESolver<M>::SolutionType;
    using typename ODE_IVPSolver<M>::Function;
    using ODE_IVPSolver<M>::ODE_IVPSolver;

    ODEStatus solve(const Function& f, const RangeType& ic) const noexcept override;
};

#include "ExplicitMidpoint.tpp"
#endif