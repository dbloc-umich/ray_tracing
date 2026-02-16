#include "FVCrankNicolson.h"
#include <iostream>

FVCrankNicolson::FVCrankNicolson(SolverPointer solver):
    _nlSolver(std::move(solver))
{
    assert(_nlSolver);
}

IVPStatus FVCrankNicolson::integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept{
    Eigen::VectorXd k1 = f(t, ic);

    auto k2Func = [&](const Eigen::VectorXd& k2) -> Eigen::VectorXd { return k2 - f(t, 0.5*dt*(k1+k2)); };
    Eigen::VectorXd k2(k1);
    _nlSolver->setFunction(k2Func);
    auto status = _nlSolver->solve(k2);
    switch (status){
        case (NLStatus::Success):
            ic += 0.5*dt*(k1+k2);
            break;
        default:
            return IVPStatus::FailureToSolve;
    }

    if (Eigen::isfinite(ic.array()).all()) return IVPStatus::Success;
    return IVPStatus::FailureToSolve;
}