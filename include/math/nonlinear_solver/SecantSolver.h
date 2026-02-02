#ifndef SECANT_SOLVER_H
#define SECANT_SOLVER_H

#include "NonlinearSolver.h"

class SecantSolver : public NonlinearSolver<1, 1>{
    public:
    using typename NonlinearSolver<1, 1>::DomainType;
    using typename NonlinearSolver<1, 1>::RangeType;
    using typename NonlinearSolver<1, 1>::Function;

    explicit SecantSolver(const Function& func, double x0, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=30);
    
    NLStatus solve(DomainType& x1) const noexcept override{ return solve(_x0, x1); }
    NLStatus solve(DomainType& x0, DomainType& x1) const noexcept;

    protected:
    mutable double _x0;
};

#endif