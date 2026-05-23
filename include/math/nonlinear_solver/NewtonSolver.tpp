#include "NewtonSolver.h"
#include "Derivative.h"
#include <iostream>

template<int N, int M>
NewtonSolver<N, M>::NewtonSolver(const Function& func, const DFunction& dfunc,
                                 double ftol, double xtol, std::size_t maxIter):
    NonlinearSolver<N, M>(func, ftol, xtol, maxIter),
    _df(dfunc)
{
    if (!_df){ // Derivative not provided, use numerical differentiation
        if constexpr(N == 1 && M == 1) _df = [this](const DomainType& x){ return df(this->_f, x); };
        else _df = [this](const DomainType& x){ return Jacobian(this->_f, x); };
    }
}

template<int N, int M>
NLStatus NewtonSolver<N, M>::solve(DomainType& x) const noexcept{
    if (!this->_f) return NLStatus::MissingFunction;

    constexpr double c = 1e-4;
    RangeType fx = this->_f(x);
    DerivativeType Jx = _df(x);

    if constexpr(N == 1){
        for (std::size_t iter = 0; iter < this->_maxIter; iter++){
            if (Jx == 0.0) return NLStatus::SingularityError;
            // std::cout << "Iter " << iter << ":" << std::endl;
            // std::cout << "x = " << x << std::endl;
            // std::cout << "fx = " << fx << std::endl;
            // std::cout << "Jx = " << Jx << std::endl;

            DomainType dx = -fx/Jx;
            double alpha; // maximum step size before hitting a bound
            if (dx > 0.0){
                if (!std::isfinite(this->_U)) alpha = 1.0;
                else{
                    double a = this->_closedU ? 1.0 : 0.99; // whether it's allowed to land on the boundary
                    alpha = std::min(1.0, a*(this->_U-x)/dx);
                }
            } else{
                if (!std::isfinite(this->_L)) alpha = 1.0;
                else{
                    double a = this->_closedL ? 1.0 : 0.99; // whether it's allowed to land on the boundary
                    alpha = std::min(1.0, a*(this->_L-x)/dx);
                }
            }

            double fx1 = this->_f(x+alpha*dx);
            while (std::abs(fx1) > std::abs(fx)*(1.0 - c*alpha) && alpha > this->_alphaMin){
                alpha = std::max(this->_alphaMin, 0.5*alpha);
                fx1 = this->_f(x+alpha*dx);
            }
            dx *= alpha;
            x += dx;
            // std::cout << "dx = " << alpha*dx << std::endl << std::endl;
            if (std::isnan(fx1)) return NLStatus::FailureToEvaluate;
            fx = fx1;

            if (this->inputConverged(x, dx) || this->outputConverged(fx)) return NLStatus::Success;
            Jx = _df(x);
        }
        return NLStatus::NoConvergence;
    }

    else{
        DerivativeType Jx = _df(x); 
        if constexpr (N == Eigen::Dynamic || M == Eigen::Dynamic){
            if (fx.size() < x.size() || fx.size() < Jx.rows() || fx.size() != Jx.cols()) return NLStatus::InvalidArgument;
            if (this->_L.size() == 0){
                this->_L = DomainType::Constant(x.size(), -mconst::infty);
                this->_closedL = BooleanType::Constant(x.size(), false);
            }
            if (this->_U.size() == 0){
                this->_U = DomainType::Constant(x.size(), mconst::infty);
                this->_closedU = BooleanType::Constant(x.size(), false);            
            }
        }

        for (std::size_t iter = 0; iter < this->_maxIter; iter++){
            Eigen::ColPivHouseholderQR<DerivativeType> qr(Jx);

            // std::cout << "Iter " << iter << ":" << std::endl;
            // std::cout << "x = " << x.transpose() << std::endl;
            // std::cout << "fx = " << fx.transpose() << std::endl;
            // std::cout << "Jx = " << std::endl << Jx << std::endl;
            if (qr.rank() < fx.size()) return NLStatus::SingularityError;

            DomainType dx = -qr.solve(fx);
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
            dx *= alpha;
            x += dx;
            // std::cout << "dx = " << dx.transpose() << std::endl << std::endl;
            if (fx1.array().isNaN().any()) return NLStatus::FailureToEvaluate;
            fx = fx1;

            if (this->inputConverged(x, dx) || this->outputConverged(fx)) return NLStatus::Success;
            Jx = _df(x);
        }  
        return NLStatus::NoConvergence;
    }    
}

template<int N, int M>
void NewtonSolver<N, M>::setFunction(const Function& f) const noexcept{
    this->_f = f;
    if constexpr (N == 1) _df = [this](const DomainType& x){ return df(this->_f, x); };
    else _df = [this](const DomainType& x){ return Jacobian(this->_f, x); };
}

template<int N, int M>
void NewtonSolver<N, M>::setFunction(const Function& f, const DFunction& df) const noexcept{
    this->_f = f;
    _df = df;
}