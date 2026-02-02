#ifndef BROYDEN_SOLVER_H
#define BROYDEN_SOLVER_H

#include "NonlinearSolver.h"

template<int N>
class BroydenSolver : public NonlinearSolver<N, N>{
    static_assert(N == Eigen::Dynamic || N > 1, "For multi-dimensional problems only. Use secant method for 1-D problems.");
    public:
    using typename NonlinearSolver<N, N>::DomainType;
    using typename NonlinearSolver<N, N>::Function;
    using DerivativeType = Eigen::Matrix<double, N, N>;
    using DFunction = std::function<DerivativeType(const DomainType&)>;

    explicit BroydenSolver(const Function& func=nullptr, const DFunction& dfuncinv=nullptr,
                           double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=20):
        NonlinearSolver<N, N>(func, ftol, xtol, maxIter){}

    NLStatus solve(DomainType& x) const noexcept override;
    NLStatus solve(DomainType& x, DerivativeType& Jinv) const noexcept;

    protected:
    mutable DFunction _dfinv; // approximation to the inverse Jacobian
};

#include "BroydenSolver.tpp"
#endif