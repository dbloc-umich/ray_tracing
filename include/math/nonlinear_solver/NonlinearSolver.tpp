#include "NonlinearSolver.h"

template<int N, int M>
NonlinearSolver<N, M>::NonlinearSolver(const Function& func, double ftol, double xtol, std::size_t maxIter):
    _f(func),
    _ftol(ftol),
    _xtol(xtol),
    _maxIter(maxIter)    
{}

template<int N, int M>
bool NonlinearSolver<N, M>::outputConverged(const RangeType& fx) const noexcept{
    if constexpr(M == 1) return std::abs(fx) <= _ftol;
    else return fx.squaredNorm() <= _ftol*_ftol;
}

template<int N, int M>
bool NonlinearSolver<N, M>::inputConverged(const DomainType& x, const DomainType& dx) const noexcept{
    if constexpr(N == 1){
        if (x == 0.0) return std::abs(dx) <= _xtol;
        return std::abs(dx/x) <= _xtol;
    }
    else{
        if (x.squaredNorm() == 0.0) return dx.squaredNorm() <= _xtol*_xtol;
        Eigen::Matrix<double, N, 1> delta(dx);
        for (Eigen::Index i = 0; i < x.size(); i++) delta[i] /= (x[i] == 0.0 ? 1 : x[i]);
        return delta.squaredNorm() <= _xtol*_xtol;
    }
}