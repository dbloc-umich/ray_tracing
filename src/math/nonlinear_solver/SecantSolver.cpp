#include "SecantSolver.h"

SecantSolver::SecantSolver(const Function& func, double x0, double ftol, double xtol, std::size_t maxIter):
    NonlinearSolver<1, 1>(func, ftol, xtol, maxIter),
    _x0(x0) // first initial guess
{};

NLStatus SecantSolver::solve(DomainType& x0, DomainType& x1) const noexcept{
    if (x0 == x1) return NLStatus::SingularityError;
    double f0 = this->_f(x0);
    double f1 = this->_f(x1);

    for (std::size_t iter = 0; iter < this->_maxIter; iter++){
        if (std::abs(f0-f1) < _ftol) return NLStatus::SingularityError;

        double tempx = x0;
        double tempf = f0;
        double dx = f0 * (x0-x1)/(f0-f1);
        x0 -= dx;
        f0 = this->_f(x0);
        if (this->inputConverged(x0, dx) || this->outputConverged(f0)) return NLStatus::Success;

        x1 = tempx;
        f1 = tempf;
    }
    return NLStatus::NoConvergence;
}