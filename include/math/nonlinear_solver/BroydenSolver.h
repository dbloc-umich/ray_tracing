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

    NLStatus solve(DomainType& x) const noexcept override{
        DerivativeType Jinv;
        if (!_dfinv) Jinv = _dfinv(x);
        else{
            if constexpr (N == 1) Jinv = DerivativeType::Identity(x.size(), x.size());
            else Jinv = DerivativeType::Identity();
        }
        return solve(x, Jinv);
    }

    NLStatus solve(DomainType& x, DerivativeType& Jinv) const noexcept{
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

    protected:
    mutable DFunction _dfinv; // approximation to the inverse Jacobian
};

#endif