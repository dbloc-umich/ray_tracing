/**
 * Solves the system of differential equations u'(t) = f(t, u), with initial condition u(0)
*/
#ifndef FV_TIME_INTEGRATOR_H
#define FV_TIME_INTEGRATOR_H

#include <functional>
#include "Eigen/Dense"

enum class IVPStatus{ Success, InvalidArgument, FailureToSolve };
class FVTimeIntegrator{
    public:
    using Function = std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>;
    virtual ~FVTimeIntegrator() = default;
    virtual IVPStatus integrate(const Function& f, Eigen::VectorXd& ic, double t, double dt) const noexcept = 0;
};

#endif