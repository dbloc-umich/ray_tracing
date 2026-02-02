#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include "Eigen/Dense"
#include <functional>
#include <type_traits>
#include <utility>

enum class ODEStatus{ Success, InvalidArgument, FailureToSolve };

template<int M>
class ODESolver{
    static_assert(M == Eigen::Dynamic || M >= 1);
    public:
    using RangeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, 1>>;
    using SolutionType = Eigen::Array<double, M, Eigen::Dynamic>;

    ODESolver() = default;
    ODESolver(double t0, double t1, double dt){ makeGrid(t0, t1, dt); }
    ODESolver(const Eigen::Array2d& t){ makeGrid(t); }
    virtual ~ODESolver() = default;

    void makeGrid(double t0, double t1, double dt) const noexcept;
    void makeGrid(const Eigen::ArrayXd& t) const noexcept;

    // Move _t and _u out of the solver object
    decltype(auto) getSolution() const noexcept{ return std::make_pair<Eigen::ArrayXd, SolutionType>(std::move(_t), std::move(_u)); } 

    protected:
    mutable Eigen::ArrayXd _t;
    mutable SolutionType _u;
};

#include "ODESolver.tpp"
#endif