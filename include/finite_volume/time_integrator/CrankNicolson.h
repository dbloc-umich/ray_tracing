#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "TimeIntegrator.h"
#include "NewtonSolver.h"

class CrankNicolson: public TimeIntegrator{
    public:
    using SolverPointer = std::unique_ptr<NonlinearSolver<Eigen::Dynamic>>;
    CrankNicolson(SolverPointer solver = std::make_unique<NewtonSolver<Eigen::Dynamic>>());   
    IVPStatus integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept override;

    protected:
    SolverPointer _nlSolver;   
};

#endif