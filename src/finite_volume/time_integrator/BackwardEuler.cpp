#include "BackwardEuler.h"

BackwardEuler::BackwardEuler(SolverPointer solver):
    _nlSolver(std::move(solver))
{
    assert(_nlSolver);
}

IVPStatus BackwardEuler::integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept{
    auto kFunc = [&](const Eigen::VectorXd& y) -> Eigen::VectorXd { return y - ic - dt*f(t,y); };
    
    Eigen::VectorXd y(ic);
    _nlSolver->setFunction(kFunc);
    auto status = _nlSolver->solve(y);
    switch (status){
        case (NLStatus::Success):
            ic += dt*y;
            break;
        default:
            return IVPStatus::FailureToSolve;
    }

    if (Eigen::isfinite(ic.array()).all()) return IVPStatus::Success;
    return IVPStatus::FailureToSolve;
}