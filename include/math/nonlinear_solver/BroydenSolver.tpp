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
    RangeType f;
    try{
        f = this->_f(x);
    } catch (...){
        return NLStatus::FailureToEvaluate;
    }

    if constexpr (N == Eigen::Dynamic){
        if (x.size() == 0 || f.size() != x.size () || f.size() != Jinv.rows() || f.size() != Jinv.cols())
            return NLStatus::InvalidArgument;
    }
    
    // double norm = this->outputNorm(f);
    // std::cout << "Initial guess: " << std::endl;
    // std::cout << "x = " << x.transpose() << ", f(x) = " << f.transpose() << std::endl;
    // std::cout << "Jinv = " << std::endl << Jinv << std::endl;

    for (std::size_t iter = 0; iter < this->_maxIter; iter++){
        DomainType dx = -Jinv*f;
        
        double alpha = 1.0;
        RangeType f1;
        try{
            f1 = this->_f(x + alpha*dx);
        } catch (const std::exception& ex){
            std::cerr << ex.what() << std::endl;
            return NLStatus::FailureToEvaluate;
        }
        // double normF1 = this->outputNorm(f1);
        // while (f1.array().isNaN().any()){
        //     std::cout << "alpha = " << alpha << std::endl;
        //     alpha *= 0.5;
        //     f1 = this->_f(x + alpha*dx);
        // }
        dx *= alpha;
        x += dx;
        // norm = normF1;
        RangeType df = f1 - f;

        // std::cout << "Iter " << iter << ":" << std::endl;
        // std::cout << "dx = " << dx.transpose() << std::endl;
        // std::cout << "x = " << x.transpose() << std::endl;
        // std::cout << "f(x) = " << f1.transpose() << std::endl;
        if (this->inputConverged(x, dx) || this->outputConverged(f1)) return NLStatus::Success;
        
        // Update the inverse Jacobian and the function value
        DomainType Jinvdf = Jinv*df;
        double denom = dx.dot(Jinvdf);
        if (denom == 0.0) return NLStatus::SingularityError;
        Eigen::Matrix<double, N, N> updateMatrix = (dx - Jinvdf)*dx.transpose() / denom;
        for (Eigen::Index i = 0; i < dx.size(); i++) updateMatrix(i,i) += 1; 
        Jinv = updateMatrix * Jinv;
        f = f1;
        // std::cout << "Jinv = " << std::endl << Jinv << std::endl << std::endl;     
    }
    return NLStatus::NoConvergence;
}