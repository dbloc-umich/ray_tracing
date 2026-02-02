#ifndef PROJECTED_NEWTON_H
#define PROJECTED_NEWTON_H

#include "ConstrainedOptimizer.h"
#include "Eigen/Cholesky"

template<int N>
class ProjectedNewton : public ConstrainedOptimizer<N>{
    public:
    using typename ConstrainedOptimizer<N>::DomainType;
    using typename ConstrainedOptimizer<N>::Function;
    using GFunction = std::function<DomainType(const DomainType&)>;
    using HessianType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, N>>;
    using HFunction = std::function<HessianType(const DomainType&)>;

    explicit ProjectedNewton(const Function& func, const DomainType& lower, const DomainType& upper,
                             const GFunction& gfunc=nullptr, const HFunction& Hfunc=nullptr,
                             double ftol=1.0e-6, double gtol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=30);

    COStatus minimize(DomainType& x) const noexcept override;
    void setFunction(const Function& f) const noexcept override;
    void setFunction(const Function& f, const GFunction& g, const HFunction& H=nullptr) const noexcept;
    void setGTol(double gtol) noexcept{ _gtol = gtol; }

    protected:
    mutable DomainType _L; // lower bound
    mutable DomainType _U; // upper bound
    mutable GFunction _g;
    mutable HFunction _H;
    double _gtol; // Tolerance in gradient

    bool gradientConverged(const DomainType& gx) const noexcept;
};

#include "ProjectedNewton.tpp"
#endif