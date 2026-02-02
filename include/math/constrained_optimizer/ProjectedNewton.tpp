#include "ProjectedNewton.h"
#include "Derivative.h"
#include "Eigen/Cholesky"

template<int N>
ProjectedNewton<N>::ProjectedNewton(const Function& func, const DomainType& lower, const DomainType& upper,
                                    const GFunction& gfunc, const HFunction& Hfunc,
                                    double ftol, double gtol, double xtol, std::size_t maxIter):
    ConstrainedOptimizer<N>(func, ftol, xtol, maxIter),
    _L(lower),
    _U(upper),
    _g(gfunc),
    _H(Hfunc),
    _gtol(gtol)
{
    if (!gfunc){
        if constexpr(N == 1) _g = [this](const DomainType& x){ return df(this->_f, x); };
        else _g = [this](const DomainType& x){ return grad(this->_f, x); };            
    }
    if (!Hfunc){
        if constexpr(N == 1) _H = [this](const DomainType& x){ return df2(this->_f, x); };
        else _H = [this](const DomainType& x){ return Hessian(this->_f, x); };
    }
}

template<int N>
COStatus ProjectedNewton<N>::minimize(DomainType& x) const noexcept{
    if (!this->_f) return COStatus::MissingFunction;

    constexpr double c = 1e-4;
    double fx = this->_f(x);
    DomainType gx = _g(x);
    HessianType Hx = _H(x);

    if constexpr(N == 1){
        if (this->gradientConverged(gx) && Hx < 0.0) return COStatus::Success;
        for (std::size_t iter = 0; iter < this->_maxIter; iter++){
            // Clamp the step size and advance
            DomainType dx = (Hx > 0.0) ? -gx/Hx : -gx;
            double alpha = dx > 0 ? std::min(1.0, (_U-x)/dx) : std::min(1.0, (_L-x)/dx); // maximum step size before hitting a bound
            double fx1 = this->_f(x+alpha*dx);
            while (fx1 > fx + c*alpha*gx*dx){
                alpha *= 0.5;
                fx1 = this->_f(x+alpha*dx);
            }
            x += alpha*dx;
            gx = _g(x);
    
            // Convergence check
            if (this->inputConverged(x, dx) || this->gradientConverged(gx) || this->outputConverged(fx, fx1)) return COStatus::Success;
            fx = fx1;
            Hx = _H(x);
        }
        return COStatus::NoConvergence;
    }

    else{
        if constexpr (N == Eigen::Dynamic){
            if (gx.size() != Hx.rows() || gx.size() != Hx.cols()) return COStatus::InvalidArgument;
        }

        for (std::size_t iter = 0; iter < this->_maxIter; iter++){
            // Check if the Hessian is positive definite
            Eigen::LLT<Eigen::Ref<HessianType>> llt(Hx);
            double lambda = 1.0e-8; // Levenberg–Marquardt damping constant
            while (llt.info() != Eigen::Success){
                Hx = llt.matrixL();
                Hx = Hx * Hx.transpose();
                for (Eigen::Index i = 0; i < gx.size(); i++) Hx(i,i) += lambda;
                lambda *= 10;
                llt.compute(Hx);
            } 

            // Clamp the step size and advance
            DomainType dx = -llt.solve(gx);
            double alpha = 1.0; // maximum step size before hitting a bound
            for (Eigen::Index i = 0; i < gx.size(); i++){
                if (dx[i] == 0.0) continue;
                alpha = dx[i] > 0 ? std::min(alpha, (_U[i]-x[i])/dx[i]) : std::min(alpha, (_L[i]-x[i])/dx[i]);
            }
            double fx1 = this->_f(x+alpha*dx);
            while (fx1 > fx + c*alpha*gx.dot(dx)){
                alpha *= 0.5;
                fx1 = this->_f(x+alpha*dx);
            }
            x += alpha*dx;
            gx = _g(x);
    
            // Convergence check
            if (this->inputConverged(x, dx) || gradientConverged(gx) || this->outputConverged(fx, fx1)) return COStatus::Success;
            fx = fx1;
            Hx = _H(x);
        }
        return COStatus::NoConvergence;
    }
}

template<int N>
void ProjectedNewton<N>::setFunction(const Function& f) const noexcept{
    ConstrainedOptimizer<N>::setFunction(f);
    if constexpr (N == 1){
        _g = [this](const DomainType& x){ return df(this->_f, x); };
        _H = [this](const DomainType& x){ return df2(this->_f, x); };
    } else{
        _g = [this](const DomainType& x){ return grad(this->_f, x); };
        _H = [this](const DomainType& x){ return Hessian(this->_f, x); };
    }
}

template<int N>
void ProjectedNewton<N>::setFunction(const Function& f, const GFunction& g, const HFunction& H) const noexcept{
    ConstrainedOptimizer<N>::setFunction(f);
    if constexpr (N == 1){
        _g = [this](const DomainType& x){ return df(this->_f, x); };
        _H = [this](const DomainType& x){ return df2(this->_f, x); };
    } else{
        _g = [this](const DomainType& x){ return grad(this->_f, x); };
        _H = [this](const DomainType& x){ return Hessian(this->_f, x); };
    }
}

template<int N>
bool ProjectedNewton<N>::gradientConverged(const DomainType& gx) const noexcept{
    if constexpr(N == 1) return std::abs(gx) <= _gtol;
    else return gx.squaredNorm() <= _gtol*_gtol; 
}