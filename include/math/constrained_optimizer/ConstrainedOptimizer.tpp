#include "ConstrainedOptimizer.h"

template<int N>
ConstrainedOptimizer<N>::ConstrainedOptimizer(const Function& func, double ftol, double xtol, std::size_t maxIter):
    _f(func),
    _ftol(ftol),
    _xtol(xtol),
    _maxIter(maxIter)
{}

template<int N>
bool ConstrainedOptimizer<N>::outputConverged(double fx0, double fx1) const noexcept{
    if (fx0 == 0.0) return std::abs(fx1) <= _ftol;
    return std::abs(1.0 - fx1/fx0) <= _ftol;
}

template<int N>
bool ConstrainedOptimizer<N>::inputConverged(const DomainType& x, const DomainType& dx) const noexcept{
    if constexpr(N == 1){
        if (x == 0) return std::abs(dx) <= _xtol;
        return std::abs(dx/x) <= _xtol;
    }
    else{
        if (x.squaredNorm() == 0.0) return dx.squaredNorm() <= _xtol*_xtol;
        Eigen::Matrix<double, N, 1> delta(dx);
        for (Eigen::Index i = 0; i < x.size(); i++) delta[i] /= (x[i] == 0.0 ? 1 : x[i]);
        return delta.squaredNorm() <= _xtol*_xtol;
    }
}
