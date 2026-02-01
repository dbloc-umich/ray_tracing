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

    void makeGrid(double t0, double t1, double dt) const noexcept{
        assert(dt != 0 && t1 != t0);
        assert((t1 > t0 && dt > 0) || (t1 < t0 && dt < 0));
        Eigen::Index steps = ceil((t1-t0)/dt);
        _t.resize(steps+1);
        _t[0] = t0;
        _t[steps] = t1;
        for (Eigen::Index i = 1; i < steps; i++) _t[i] = _t[i-1] + dt;
        _u.resize(M == Eigen::Dynamic ? 1 : M, steps+1);
    }

    void makeGrid(const Eigen::ArrayXd& t){
        assert(std::is_sorted(t.cbegin(), t.cend(), [](double a, double b){ return a < b; } ));
        _t = t;
        _u.resize(M == Eigen::Dynamic ? 1 : M, _t.size());
    }

    decltype(auto) getSolution() const noexcept{
        // Move _t and _u out of the solver object
        return std::make_pair<Eigen::ArrayXd, SolutionType>(std::move(_t), std::move(_u));
    } 

    protected:
    mutable Eigen::ArrayXd _t;
    mutable SolutionType _u;
};

#endif