#include "ForwardEuler.h"
#include <iostream>

IVPStatus ForwardEuler::integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const{
    try{
        u0 += dt*f(t, u0);
        if (Eigen::isfinite(u0.array()).all()) return IVPStatus::Success;
        return IVPStatus::FailureToEvaluate;
    } catch(const std::runtime_error& ex){
        std::cerr << ex.what() << std::endl;
        return IVPStatus::FailureToEvaluate;
    }
}