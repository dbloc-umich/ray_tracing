#include "BroydenSolver.h"

template<int N>
NLStatus BroydenSolver<N>::solve(DomainType& x) const noexcept{
    DerivativeType Jinv;
    if (!_dfinv) Jinv = _dfinv(x);
    else{
        if constexpr (N == 1) Jinv = DerivativeType::Identity(x.size(), x.size());
        else Jinv = DerivativeType::Identity();
    }
    return solve(x, Jinv);
}

template<int N>
NLStatus BroydenSolver<N>::solve(DomainType& x, DerivativeType& Jinv) const noexcept{
    if (!this->_f) return NLStatus::MissingFunction;
    DomainType f = this->_f(x);
    if (N == Eigen::Dynamic){
        if (x.size() == 1 || f.size() != x.size () || f.size() != Jinv.rows() || f.size() != Jinv.cols())
            return NLStatus::InvalidArgument;
    }

    for (std::size_t iter = 0; iter < this->_maxIter; iter++){
        DomainType dx = Jinv*f;
        x -= dx;
        DomainType f1 = this->_f(x);
        DomainType df = f1 - f;
        if (this->inputConverged(x, dx) || this->outputConverged(f1)) return NLStatus::Success;
        
        // Update the inverse Jacobian and the function value
        DomainType Jinvdf = Jinv*df;
        double denom = dx.dot(Jinvdf);
        if (denom == 0.0) return NLStatus::SingularityError;
        Eigen::Matrix<double, N, N> updateMatrix = (dx - Jinvdf).outer(dx) / denom;
        for (Eigen::Index i = 0; i < dx.size(); i++) updateMatrix(i,i) += 1; 
        Jinv = updateMatrix * Jinv;
    }
    return NLStatus::NoConvergence;
}