#include "DIRK2.h"
#include "NonlinearSolver.h"
#include <iostream>

IVPStatus DIRK2::integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const{
    if (!_nlSolver) return IVPStatus::MissingNonlinearSolver;

    auto k1Func = [&, this](const Eigen::VectorXd& k1) -> Eigen::VectorXd { return k1 - f(t+_x*dt, u0 + dt*_x*k1); };
    _nlSolver->setFunction(k1Func);
    Eigen::VectorXd k1 = f(t, u0);
    auto status = _nlSolver->solve(k1);
    switch (status){
        case (NLStatus::Success):
            if (!Eigen::isfinite(k1.array()).all()) return IVPStatus::FailureToEvaluate;
            break;
        default:
            return IVPStatus::FailureToSolve;
    }

    auto k2Func = [&, this](const Eigen::VectorXd& k2) -> Eigen::VectorXd { return k2 - f(t+dt, u0 + dt*((1-_x)*k1 + _x*k2)); };
    _nlSolver->setFunction(k2Func);
    Eigen::VectorXd k2(k1);
    status = _nlSolver->solve(k2);
    switch (status){
        case (NLStatus::Success):
            if (!Eigen::isfinite(k2.array()).all()) return IVPStatus::FailureToEvaluate;
            u0 += dt*((1-_x)*k1 + _x*k2);
            return IVPStatus::Success;
        default:
            return IVPStatus::FailureToSolve;
    }    
}