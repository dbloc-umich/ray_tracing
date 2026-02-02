#ifndef CONSTRAINED_OPTIMIZER_H
#define CONSTRAINED_OPTIMIZER_H

#include "Eigen/Dense"
#include <functional>
#include <type_traits>

enum class COStatus{ Success, MissingFunction, InvalidArgument, NoConvergence };

template<int N>
class ConstrainedOptimizer{
    static_assert(N == Eigen::Dynamic || N >= 1);
    public:
    using DomainType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, 1>>;
    using Function = std::function<double(const DomainType&)>;

    explicit ConstrainedOptimizer(const Function& func=nullptr, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=100);
    virtual ~ConstrainedOptimizer() = default;

    virtual COStatus minimize(DomainType&) const noexcept = 0;
    virtual void setFunction(const Function& f) const noexcept{ _f = f; }

    void setFTol(double ftol) noexcept{ _ftol = ftol; }
    void setXTol(double xtol) noexcept{ _xtol = xtol; }
    void setMaxIter(std::size_t maxIter) noexcept{ _maxIter = maxIter; }

    protected:
    mutable Function _f;
    double _ftol; // Tolerance in output (objective function) step size
    double _xtol; // Tolerance in input step size
    std::size_t _maxIter;

    bool outputConverged(double fx0, double fx1) const noexcept;
    bool inputConverged(const DomainType& x, const DomainType& dx) const noexcept;
};

#include "ConstrainedOptimizer.tpp"
#endif