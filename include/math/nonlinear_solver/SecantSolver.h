#ifndef SECANT_SOLVER_H
#define SECANT_SOLVER_H

#include "NonlinearSolver.h"

class SecantSolver : public NonlinearSolver<1, 1>{
    public:
    using typename NonlinearSolver<1, 1>::DomainType;
    using typename NonlinearSolver<1, 1>::RangeType;
    using typename NonlinearSolver<1, 1>::Function;

    explicit SecantSolver(const Function& func, double x0, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=30);
    
    // x1 is the better guess out of the two guesses
    NLStatus solve(DomainType& x1) const noexcept override{ return solve(x1, _x0); }
    NLStatus solve(DomainType& x1, DomainType& x0) const noexcept;

    protected:
    mutable double _x0;
};

#endif