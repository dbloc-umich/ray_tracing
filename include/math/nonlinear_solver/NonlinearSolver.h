#ifndef NONLINEAR_SOLVER_H
#define NONLINEAR_SOLVER_H

#include <Eigen/Dense>
#include <functional>
#include <type_traits>

enum class NLStatus{ Success, MissingFunction, InvalidArgument, SingularityError, NoConvergence, FailureToEvaluate };

template<int N, int M = N>
class NonlinearSolver{
    static_assert(N == Eigen::Dynamic || N >= 1);
    static_assert(M == Eigen::Dynamic || M >= N);
    public:
    using DomainType = std::conditional_t<N == 1, double, Eigen::Vector<double, N>>;
    using RangeType = std::conditional_t<M == 1, double, Eigen::Vector<double, M>>;
    using BooleanType = std::conditional_t<N == 1, bool, Eigen::Vector<bool, N>>;
    using Function = std::function<RangeType(const DomainType&)>;

    explicit NonlinearSolver(const Function& func=nullptr, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=100);
    virtual ~NonlinearSolver() = default;

    virtual NLStatus solve(DomainType&) const noexcept = 0;
    virtual void setFunction(const Function& f) const noexcept{ _f = f; }

    template<int NN = N, typename = std::enable_if_t<NN != 1>>
    double lowerBound(Eigen::Index i) const noexcept{ return _L[i]; } 
    const DomainType& lowerBound() const noexcept{ return _L; }

    template<int NN = N, typename = std::enable_if_t<NN != 1>>
    void setLowerBound(Eigen::Index i, double L, bool isClosed = false);
    void setLowerBound(const DomainType& _L, const BooleanType& isClosed = defaultBooleanType());

    template<int NN = N, typename = std::enable_if_t<NN != 1>>
    double upperBound(Eigen::Index i) const noexcept{ return _U[i]; }
    const DomainType& upperBound() const noexcept{ return _U; }

    template<int NN = N, typename = std::enable_if_t<NN != 1>>
    void setUpperBound(Eigen::Index i, double U, bool isClosed = false);
    void setUpperBound(const DomainType& U, const BooleanType& isClosed = defaultBooleanType());

    void setAlphaMin(double alphaMin);
    void setFTol(double ftol);
    void setXTol(double xtol);
    void setMaxIter(std::size_t maxIter);

    protected:
    mutable Function _f;

    // Line search parameters
    mutable DomainType _L;
    mutable BooleanType _closedL;
    mutable DomainType _U;
    mutable BooleanType _closedU;
    double _alphaMin;

    // Convergence parameters    
    double _ftol; // Tolerance in output
    double _xtol; // Tolerance in input step size
    std::size_t _maxIter;

    bool outputConverged(const RangeType& fx) const noexcept;
    bool inputConverged(const DomainType& x, const DomainType& dx) const noexcept;
    static BooleanType defaultBooleanType();
};

#include "NonlinearSolver.tpp"
#endif
