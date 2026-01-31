#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "ODE_IVPSolver.h"
#include "NonlinearSolver.h"

template<int M, typename... Args>
class CrankNicolson: public ODE_IVPSolver<M>{
    public:
    using typename ODESolver<M>::RangeType;
    using typename ODESolver<M>::SolutionType;
    using typename ODE_IVPSolver<M>::Function;
    using SolverPointer = std::unique_ptr<NonlinearSolver<M, M, Args...>>;
    
    CrankNicolson(double t0, double t1, double dt, SolverPointer solver=nullptr):
        ODE_IVPSolver(t0, t1, dt),
        _nlSolver(std::move(solver))
    {}
    CrankNicolson(double t0, double t1, std::size_t steps, SolverPointer solver=nullptr):
        ODE_IVPSolver(t0, t1, steps),
        _nlSolver(std::move(solver))
    {}
    CrankNicolson(const Eigen::ArrayXd& t, SolverPointer solver=nullptr):
        ODE_IVPSolver(t),
        _nlSolver(std::move(solver))
    {}    

    protected:
    SolverPointer _nlSolver;   
};

#endif