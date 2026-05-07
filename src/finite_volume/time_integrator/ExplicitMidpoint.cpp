#include "ExplicitMidpoint.h"
#include <iostream>

IVPStatus ExplicitMidpoint::integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept{
    try{
        // std::cout << "u0 = " << ic.transpose() << std::endl;
        Eigen::VectorXd k1 = f(t, ic);
        // std::cout << "k1 = " << k1.transpose() << std::endl;
        // std::cout << "u1/2 = " << (ic+0.5*dt*k1).transpose() << std::endl;
        Eigen::VectorXd k2 = f(t+0.5*dt, ic+0.5*dt*k1);
        // std::cout << "k2 = " << k2.transpose() << std::endl;
        ic += dt*k2;
        // std::cout << "u1 = " << ic.transpose() << std::endl;
    } catch(const std::runtime_error& ex){
        std::cerr << ex.what() << std::endl;
        return IVPStatus::FailureToSolve;
    }

    if (Eigen::isfinite(ic.array()).all()) return IVPStatus::Success;
    return IVPStatus::FailureToSolve;
}