#ifndef FV_CRANK_NICOLSON_H
#define FV_CRANK_NICOLSON_H

#include "FVTimeIntegrator.h"
#include "NewtonSolver.h"

class FVCrankNicolson: public FVTimeIntegrator{
    public:
    using SolverPointer = std::unique_ptr<NonlinearSolver<Eigen::Dynamic>>;
    FVCrankNicolson(SolverPointer solver = std::make_unique<NewtonSolver<Eigen::Dynamic>>());   
    IVPStatus integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept override;

    protected:
    SolverPointer _nlSolver;   
};

#endif