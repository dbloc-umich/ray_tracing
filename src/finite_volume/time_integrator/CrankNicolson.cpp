#include "CrankNicolson.h"
#include "NonlinearSolver.h"
#include <iostream>

IVPStatus CrankNicolson::integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const{
    if (!_nlSolver) return IVPStatus::MissingNonlinearSolver;

    auto uFunc = [&](const Eigen::VectorXd& u) -> Eigen::VectorXd { return u - u0 - 0.5*dt*(f(t,u0) + f(t+dt,u)); };
    _nlSolver->setFunction(uFunc);
    Eigen::VectorXd u(u0);
    
    auto status = _nlSolver->solve(u);
    switch (status){
        case (NLStatus::Success):
            if (Eigen::isfinite(u.array()).all()){
                u0 = std::move(u);
                return IVPStatus::Success;
            }
            return IVPStatus::FailureToEvaluate;
        default:
            return IVPStatus::FailureToSolve;
    }
}