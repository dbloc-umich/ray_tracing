#include "BroydenSolver.h"

template<int N>
NLStatus BroydenSolver<N>::solve(DomainType& x) const noexcept{
    DerivativeType Jinv;
    if (_dfinv) Jinv = _dfinv(x);
    else{
        if constexpr (N == Eigen::Dynamic) Jinv = DerivativeType::Identity(x.size(), x.size());
        else Jinv = DerivativeType::Identity();
    }
    return solve(x, Jinv);
}

template<int N>
NLStatus BroydenSolver<N>::solve(DomainType& x, DerivativeType& Jinv) const noexcept{
    if (!this->_f) return NLStatus::MissingFunction;
    constexpr double c = 1e-4;
    RangeType fx = this->_f(x);

    if constexpr (N == Eigen::Dynamic){
        if (x.size() == 0 || fx.size() != x.size () || fx.size() != Jinv.rows() || fx.size() != Jinv.cols())
            return NLStatus::InvalidArgument;
        if (this->_L.size() == 0){
            this->_L = DomainType::Constant(x.size(), -mconst::infty);
            this->_closedL = BooleanType::Constant(x.size(), false);
        }
        if (this->_U.size() == 0){
            this->_U = DomainType::Constant(x.size(), mconst::infty);
            this->_closedU = BooleanType::Constant(x.size(), false);            
        }
    }
    
    // std::cout << "Initial guess: " << std::endl;
    // std::cout << "x = " << x.transpose() << std::endl;
    // std::cout << "f(x) = " << fx.transpose() << std::endl;
    // std::cout << "Jinv = " << std::endl << Jinv << std::endl;

    for (std::size_t iter = 0; iter < this->_maxIter; iter++){
        // std::cout << "Iter " << iter << ":" << std::endl;
        DomainType dx = -Jinv*fx;
        double alpha = 1.0; // maximum step size before hitting a bound
        for (Eigen::Index i = 0; i < x.size(); i++){
            if (dx[i] == 0.0) continue;
            else if (dx[i] > 0.0){
                if (!std::isfinite(this->_U[i])) alpha = std::min(alpha, 1.0);
                else{
                    double a = this->_closedU[i] ? 1.0 : 0.99; // whether it's allowed to land on the boundary
                    alpha = std::min(alpha, a*(this->_U[i]-x[i])/dx[i]);
                }
            } else{
                if (!std::isfinite(this->_L[i])) alpha = std::min(alpha, 1.0);
                else{
                    double a = this->_closedL[i] ? 1.0 : 0.99; // whether it's allowed to land on the boundary
                    alpha = std::min(alpha, a*(this->_L[i]-x[i])/dx[i]);
                }
            }
        }
        // std::cout << "Initial alpha = " << alpha << std::endl;
        // std::cout << "x+alpha*dx = " << (x+alpha*dx).transpose() << std::endl;
        RangeType fx1 = this->_f(x+alpha*dx);
        // std::cout << "fx1 = " << fx1.transpose() << std::endl;
        while (fx1.norm() > fx.norm()*(1.0 - c*alpha) && alpha > this->_alphaMin){
            alpha = std::max(this->_alphaMin, 0.5*alpha);
            // std::cout << "next alpha = " << alpha << std::endl;
            // std::cout << "x+alpha*dx = " << (x+alpha*dx).transpose() << std::endl;
            fx1 = this->_f(x+alpha*dx);
            // std::cout << "fx1 = " << fx1.transpose() << std::endl;
        }
        // std::cout << "Final alpha = " << alpha << std::endl;
        dx *= alpha;
        x += dx;
        if (fx1.array().isNaN().any()) return NLStatus::FailureToEvaluate;
        RangeType df = fx1 - fx;

        // std::cout << "dx = " << dx.transpose() << " with alpha = " << alpha << std::endl;
        // std::cout << "x = " << x.transpose() << std::endl;
        // std::cout << "f(x) = " << fx1.transpose() << std::endl;
        if (this->inputConverged(x, dx) || this->outputConverged(fx1)) return NLStatus::Success;
        
        // Update the inverse Jacobian and the function value
        DomainType Jinvdf = Jinv*df;
        double denom = dx.dot(Jinvdf);
        if (denom == 0.0) return NLStatus::SingularityError;
        Eigen::Matrix<double, N, N> updateMatrix = (dx - Jinvdf)*dx.transpose() / denom;
        for (Eigen::Index i = 0; i < dx.size(); i++) updateMatrix(i,i) += 1; 
        Jinv = updateMatrix * Jinv;
        fx = fx1;
        // std::cout << "Jinv = " << std::endl << Jinv << std::endl << std::endl;     
    }
    // std::cout << std::endl;
    return NLStatus::NoConvergence;
}