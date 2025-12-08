#ifndef NONLINEAR_SOLVER_H
#define NONLINEAR_SOLVER_H

#include <Eigen/Dense>
#include <functional>
#include <type_traits>

template<int Domain, int Range = Domain,
         typename = std::enable_if_t<(Domain == -1 || Domain >= 1) && (Range == -1 || Range >= Domain)>>
class NonlinearSolver{
    public:
    enum class Status{ Success, InvalidGuess, InvalidArgument, SolverError, NoConvergence };
    using DomainType = std::conditional_t<Domain == 1, double, Eigen::Matrix<double, Domain, 1>>;
    using RangeType = std::conditional_t<Range == 1, double, Eigen::Matrix<double, Range, 1>>;
    using Function = std::function<RangeType(const DomainType&)>;

    explicit NonlinearSolver(const Function& func=nullptr, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=100):
        _f(func),
        _ftol(ftol),
        _xtol(xtol),
        _maxIter(maxIter)    
    {}
    virtual ~NonlinearSolver() = default;

    virtual Status solve(DomainType&) const noexcept = 0;
    Status solve(DomainType&& x0) const noexcept{
        DomainType x(std::move(x0));
        Status stat = solve(f, x);
        _result = std::move(x); // store the result in a member variable
        return stat;
    }
    const DomainType& result() const noexcept{ return _result; }

    void setFunction(const Function& f) const noexcept{ _f = f; }

    protected:
    mutable Function _f;
    double _ftol; // Tolerance in output
    double _xtol; // Tolerance in input between two successive iterations
    std::size_t _maxIter;
    mutable DomainType _result; // stores the result of the solver

    bool outputConverged(const RangeType& fx){
        if constexpr(Range == 1) return std::abs(fx) <= _ftol;
        return fx.squaredNorm() <= _ftol*_ftol;
    }

    bool inputConverged(const DomainType& x, const DomainType& dx){
        if constexpr(Domain == 1){
            if (x == 0) return std::abs(dx) <= _xtol;
            return std::abs(dx/x) <= _xtol;
        }
        if (x.squaredNorm() == 0.0) return dx.squaredNorm() <= _xtol*_xtol;
        Eigen::Matrix<double, D, 1> delta = dx;
        for (Eigen::Index i = 0; i < x.size(); i++) delta[i] /= (x[i] == 0.0 ? 1 : x[i]);
        return delta.squaredNorm() <= _xtol*_xtol;
    }
};
#endif