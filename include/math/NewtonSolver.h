#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include "NonlinearSolver.h"
#include "Derivative.h"

template<int Domain, int Range = Domain,
         typename = std::enable_if_t<(Domain == -1 || Domain >= 1) && (Range == -1 || Range >= Domain)>>
class NewtonSolver : NonlinearSolver<Domain, Range>{
    public:
    using typename NonlinearSolver<Domain, Range>::DomainType;
    using typename NonlinearSolver<Domain, Range>::RangeType;
    using typename NonlinearSolver<Domain, Range>::Function;
    using DerivativeType = std::conditional_t<Range == 1, double, Eigen::Matrix<double, Range, Domain>>;
    using DFunction = std::function<DerivativeType(const DomainType&)>;

    explicit NewtonSolver(const Function& func=nullptr, const DFunction& dfunc=nullptr,
                          double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=100):
        NonlinearSolver<Domain, Range>(func, ftol, xtol, maxIter),
        _df(dfunc)
    {
        if (!_df){ // Derivative not provided, use numerical differentiation
            if constexpr(Domain == 1 && Range == 1) _df = [this](const DomainType& x){ return df(_f, x); };
            _df = [this](const DomainType& x){ return Jacobian(_f, x); };
        }
    }

    Status solve(DomainType& x) const noexcept override{
        if constexpr(Domain == 1){
            for (std::size_t iter = 0; iter < _maxIter; iter++){
                DerivativeType dfx = _df(x);
                if (dfx == 0.0) return Status::SingularityError;

                DomainType dx = -f(x)/dfx;
                x += dx;
                if (inputConverged(x, dx)) return Status::Success;
            }
            return Status::NoConvergence;
        }

        bool dimensionsChecked = false;
        for (std::size_t iter = 0; iter < maxIter; iter++){
            RangeType fx = _f(x);
            DerivativeType Jx = _df(x);
            // std::cout << "iter = " << iter << ": x = " << x.transpose() << ", f(x) = " << fx.transpose() << std::endl;
            // std::cout << "Jx = " << std::endl << Jx << std::endl;
            if (!dimensionsChecked){
                // Check against dimension mismatch, for the first iteration only
                if (fx.size() < x.size() || fx.size() < Jx.rows()) return Status::InvalidArgument;
                if (fx.size() != Jx.cols()) return Status::InvalidArgument;
                dimensionsChecked = true;
            }
            Eigen::ColPivHouseholderQR<Eigen::Ref<DerivativeType> > > qr(Jx);
            if (qr.rank() < fx.size()) return Status::SingularityError;
            auto dx = qr.solve(fx);
            x -= dx;

            // std::cout << "Error = " << dx.squaredNorm() / x.squaredNorm() << std::endl << std::endl;
            if (inputConverged(x, dx)) return Status::Success;
        }
        return Status::NoConvergence;
    }

    protected:
    mutable DFunction _df;
};

#endif