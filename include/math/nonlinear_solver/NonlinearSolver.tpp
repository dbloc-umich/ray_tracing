#include "NonlinearSolver.h"
#include "Constants.h"
#include <iostream>

template<int N, int M>
NonlinearSolver<N, M>::NonlinearSolver(const Function& func, double ftol, double xtol, std::size_t maxIter):
    _f(func),
    _alphaMin(0.01),
    _ftol(ftol),
    _xtol(xtol),
    _maxIter(maxIter)    
{
    if constexpr (N == 1){
        _L = -mconst::infty;
        _U =  mconst::infty;
    } else if constexpr (N != Eigen::Dynamic){
        _L = DomainType::Constant(-mconst::infty);
        _U = DomainType::Constant( mconst::infty);
    }
}

template<int N, int M>
template<int NN, typename>
void NonlinearSolver<N, M>::setLowerBound(Eigen::Index i, double L, bool isClosed){
    if (_U[i] > L){
        _L[i] = L;
        _closedL[i] = std::isfinite(L) ? isClosed : false;
    }
    else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
} 

template<int N, int M>
void NonlinearSolver<N, M>::setLowerBound(const DomainType& L, const BooleanType& isClosed){
    if constexpr(N == 1){
        if (_U > L){
            _L = L;
            _closedL = std::isfinite(L) ? isClosed : false;
        } else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
    } else if constexpr(N == Eigen::Dynamic){
        assert(L.size() == _L.size() || _L.size() == 0);
        assert(isClosed.size() == _closedL.size() || isClosed.size() == 0 || _closedL.size() == 0);

        if (_L.size() == 0) _L = DomainType::Constant(L.size(), -mconst::infty);
        if (_closedL.size() == 0) _closedL = BooleanType::Constant(L.size(), false);
        if (_U.size() == 0) _U = DomainType::Constant(L.size(), mconst::infty);
        if (_closedU.size() == 0) _closedU = BooleanType::Constant(L.size(), false);        

        for (Eigen::Index i = 0; i < _L.size(); i++){
            if (isClosed.size() == 0) setLowerBound(i, L[i], false);
            else setLowerBound(i, L[i], isClosed[i]);
        }
    } else{
        for (Eigen::Index i = 0; i < _L.size(); i++) setLowerBound(i, L[i], isClosed[i]);
    }
}

template<int N, int M>
template<int NN, typename>
void NonlinearSolver<N, M>::setUpperBound(Eigen::Index i, double U, bool isClosed){
    if (_L[i] < U){
        _U[i] = U;
        _closedU[i] = std::isfinite(U) ? isClosed : false;
    }
    else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
} 

template<int N, int M>
void NonlinearSolver<N, M>::setUpperBound(const DomainType& U, const BooleanType& isClosed){
    if constexpr(N == 1){
        if (_L < U){
            _U = U;
            _closedU = std::isfinite(U) ? isClosed : false;
        } else throw std::invalid_argument("ERROR: Lower bound must be strictly less than upper bound.");
    } else if constexpr(N == Eigen::Dynamic){
        assert(U.size() == _U.size() || _U.size() == 0);
        assert(isClosed.size() == _closedU.size() || isClosed.size() == 0 || _closedU.size() == 0);

        if (_L.size() == 0) _L = DomainType::Constant(U.size(), -mconst::infty);
        if (_closedL.size() == 0) _closedL = BooleanType::Constant(U.size(), false);
        if (_U.size() == 0) _U = DomainType::Constant(U.size(), mconst::infty);
        if (_closedU.size() == 0) _closedU = BooleanType::Constant(U.size(), false);

        for (Eigen::Index i = 0; i < _U.size(); i++){
            if (isClosed.size() == 0) setUpperBound(i, U[i], false);
            else setUpperBound(i, U[i], isClosed[i]);
        }
    } else{
        for (Eigen::Index i = 0; i < _U.size(); i++) setUpperBound(i, U[i], isClosed[i]);
    }
}
template<int N, int M>
void NonlinearSolver<N, M>::setAlphaMin(double alphaMin){
    if (alphaMin <= 0.0 || alphaMin > 1) throw std::invalid_argument("ERROR: alphaMin must be in the interval (0, 1].");
    _alphaMin = alphaMin;
}

template<int N, int M>
void NonlinearSolver<N, M>::setFTol(double ftol){
    if (ftol <= 0.0) throw std::invalid_argument("ERROR: Positive ftol required.");
    _ftol = ftol;
}

template<int N, int M>
void NonlinearSolver<N, M>::setXTol(double xtol){
    if (xtol <= 0.0) throw std::invalid_argument("ERROR: Positive xtol required.");
    _xtol = xtol;    
}

template<int N, int M>
void NonlinearSolver<N, M>::setMaxIter(std::size_t maxIter){
    if (maxIter == 0) throw std::invalid_argument("ERROR: Positive maxIter required.");
    _maxIter = maxIter;
}

template<int N, int M>
bool NonlinearSolver<N, M>::outputConverged(const RangeType& fx) const noexcept{
    if constexpr(M == 1) return std::abs(fx) < _ftol;
    else return fx.squaredNorm() < _ftol*_ftol;
}

template<int N, int M>
bool NonlinearSolver<N, M>::inputConverged(const DomainType& x, const DomainType& dx) const noexcept{
    if constexpr(N == 1){
        if (std::abs(x) < _xtol) return std::abs(dx) < _xtol;
        return std::abs(dx/x) < _xtol;
    } else{
        if (x.squaredNorm() == 0.0) return dx.squaredNorm() < _xtol*_xtol;
        Eigen::Matrix<double, N, 1> delta(dx);
        for (Eigen::Index i = 0; i < x.size(); i++) delta[i] /= (std::fabs(x[i]) < _xtol ? 1 : x[i]);
        return delta.squaredNorm() < _xtol*_xtol;
    }
}

template<int N, int M>
typename NonlinearSolver<N, M>::BooleanType NonlinearSolver<N, M>::defaultBooleanType(){
    if constexpr(N == 1) return false;
    else if constexpr(N == Eigen::Dynamic) return BooleanType::Constant(0, false);
    else return BooleanType::Constant(false);
}