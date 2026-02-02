#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include "NonlinearSolver.h"

template<int N, int M = N>
class NewtonSolver : public NonlinearSolver<N, M>{
    public:
    using typename NonlinearSolver<N, M>::DomainType;
    using typename NonlinearSolver<N, M>::RangeType;
    using typename NonlinearSolver<N, M>::Function;
    using DerivativeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, N>>;
    using DFunction = std::function<DerivativeType(const DomainType&)>;

    explicit NewtonSolver(const Function& func=nullptr, const DFunction& dfunc=nullptr,
                          double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=20);

    NLStatus solve(DomainType& x) const noexcept override;
    void setFunction(const Function& f) const noexcept override;
    void setFunction(const Function& f, const DFunction& df) const noexcept;

    protected:
    mutable DFunction _df;
};

#include "NewtonSolver.tpp"
#endif