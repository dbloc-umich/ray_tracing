#include "SecantSolver.h"
#include <iostream>

SecantSolver::SecantSolver(const Function& func, double x0, double ftol, double xtol, std::size_t maxIter):
    NonlinearSolver<1, 1>(func, ftol, xtol, maxIter),
    _x0(x0) // first initial guess
{};

NLStatus SecantSolver::solve(DomainType& x1, DomainType& x0) const noexcept{
    if (x0 == x1) return NLStatus::SingularityError;
    constexpr double c = 1e-4;
    double f0 = this->_f(x0);
    double f1 = this->_f(x1);
    std::cout << "Initial guess: x0 = " << x0 << ", x1 = " << x1 << std::endl;

    for (std::size_t iter = 0; iter < this->_maxIter; iter++){
        if (std::abs(f0-f1) < _ftol) return NLStatus::SingularityError;

        double tempx = x1;
        double tempf = f1;

        double dx = -f1 * (x1-x0)/(f1-f0);
        double alpha = 1.0; // maximum step size before hitting a bound
        if (dx > 0.0){
            if (!std::isfinite(this->_U)) alpha = 1.0;
            else{
                double a = this->_closedU ? 1.0 : 0.99; // whether it's allowed to land on the boundary
                alpha = std::min(alpha, a*(this->_U-x1)/dx);
            }
        } else{
            if (!std::isfinite(this->_L)) alpha = 1.0;
            else{
                double a = this->_closedL ? 1.0 : 0.99; // whether it's allowed to land on the boundary
                alpha = std::min(alpha, a*(this->_L-x1)/dx);
            }
        }

        double f2 = this->_f(x1+alpha*dx);
        while (std::abs(f2) > std::abs(f1)*(1.0 - c*alpha) && alpha > this->_alphaMin){
            alpha = std::max(this->_alphaMin, 0.5*alpha);
            f2 = this->_f(x1+alpha*dx);
        }
        dx *= alpha;
        x1 += dx;
        if (std::isnan(f2)) return NLStatus::FailureToEvaluate;
        f1 = f2;
        if (this->inputConverged(x1, dx) || this->outputConverged(f1)) return NLStatus::Success;

        x0 = tempx;
        f0 = tempf;
        std::cout << "Iter " << iter << ": x0 = " << x0 << ", x1 = " << x1 << std::endl;
    }
    return NLStatus::NoConvergence;
}