/**
 * Solves the system of differential equations u'(t) = f(t, u)
*/
#ifndef ODE_IVP_SOLVER_H
#define ODE_IVP_SOLVER_H

#include "ODESolver.h"

template<int M>
class ODE_IVPSolver: public ODESolver<M>{
    public:
    using typename ODESolver<M>::RangeType;
    using typename ODESolver<M>::SolutionType;
    using Function = std::function<RangeType(double, const RangeType&)>;
    using ODESolver<M>::ODESolver;

    virtual ~ODE_IVPSolver() = default;
    virtual ODEStatus solve(const Function& f, const RangeType& ic) const noexcept = 0;
};

#endif