#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "ODE_IVPSolver.h"
#include "NewtonSolver.h"

template<int M>
class CrankNicolson: public ODE_IVPSolver<M>{
    public:
    using typename ODESolver<M>::RangeType;
    using typename ODESolver<M>::SolutionType;
    using typename ODE_IVPSolver<M>::Function;
    using SolverPointer = std::unique_ptr<NonlinearSolver<M>>;
    
    CrankNicolson(double t0, double t1, double dt, SolverPointer solver=std::make_unique<NewtonSolver<M>>());
    CrankNicolson(const Eigen::ArrayXd& t, SolverPointer solver=std::make_unique<NewtonSolver<M>>());   

    ODEStatus solve(const Function& f, const RangeType& ic) const noexcept override;

    protected:
    SolverPointer _nlSolver;   
};

#include <CrankNicolson.tpp>
#endif