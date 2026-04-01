#include "BackwardEuler.h"
#include <iostream>

BackwardEuler::BackwardEuler(SolverPointer solver):
    _nlSolver(std::move(solver))
{
    assert(_nlSolver);
}

IVPStatus BackwardEuler::integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const noexcept{
    auto kFunc = [&](const Eigen::VectorXd& y) -> Eigen::VectorXd { return y - u0 - dt*f(t,y); };
    
    Eigen::VectorXd y(u0);
    _nlSolver->setFunction(kFunc);
    auto status = _nlSolver->solve(y);
    switch (status){
        case (NLStatus::Success):
            u0 = std::move(y);
            break;
        default:
            return IVPStatus::FailureToSolve;
    }

    if (Eigen::isfinite(u0.array()).all()) return IVPStatus::Success;
    return IVPStatus::FailureToSolve;
}