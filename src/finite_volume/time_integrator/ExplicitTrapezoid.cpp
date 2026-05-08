#include "ExplicitTrapezoid.h"
#include <iostream>

IVPStatus ExplicitTrapezoid::integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const{
    try{
        Eigen::VectorXd k1 = f(t, u0);
        Eigen::VectorXd u1 = u0 + dt*k1;
        u0 += 0.5*dt*(k1 + f(t+dt, u1));
        if (Eigen::isfinite(u0.array()).all()) return IVPStatus::Success;
        return IVPStatus::FailureToEvaluate;
    } catch(const std::runtime_error& ex){
        std::cerr << ex.what() << std::endl;
        return IVPStatus::FailureToEvaluate;
    }
}