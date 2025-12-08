// #ifndef NONLINEAR_SOLVERS
// #define NONLINEAR_SOLVERS

// #include <Eigen/Dense>
// #include <functional>
// #include <type_traits>

// template<int Domain, int Range = Domain,
//          typename = std::enable_if_t<(Domain == -1 || Domain >= 1) && (Range == -1 || Range >= Domain)>>
// class NonlinearSolver{
//     public:
//     enum class Status{ Success, InvalidGuess, SolverError, NoConvergence }
//     using DomainType = std::conditional_t<Domain == 1, double, Eigen::Matrix<double, Domain, 1>>;
//     using RangeType = std::conditional_t<Range == 1, double, Eigen::Matrix<double, Range, 1>>;
//     using Function = std::function<RangeType(const DomainType&)>&;

//     NonlinearSolver(double tol, std::size_t maxIter): _tol(tol), _maxIter(maxIter) {}
//     virtual ~NonlinearSolver() = default;

//     virtual Status solve(const &) const noexcept = 0;
//     Status solve(const std::function<RangeType(const DomainType&)>& f, DomainType&& x0) const noexcept{
//         DomainType x(std::move(x0));
//         Status stat = solve(f, x);
//         _result = x; // store the result in a member variable
//         return stat;
//     }

//     DomainType& result() const noexcept{ return _result; }

//     private:
//     double _tol;
//     std::size_t _maxIter;
//     mutable DomainType _result; // stores the result of the solver
// };
// #endif

// #include "Eigen/Dense"
// #include "Derivative.h"
// #include <cmath>
// #include <type_traits>

// // Failure cases
// class NewtonSingularityError : public std::runtime_error {
// public:
//     explicit NewtonSingularityError(const std::string& msg) : std::runtime_error(msg) {}
// };


// // Thrown when Newton's method exceeds maximum allowed iterations
// class NewtonMaxIterationsError : public std::runtime_error {
// public:
//     explicit NewtonMaxIterationsError(const std::string& msg) : std::runtime_error(msg) {}
// };

// // Newton's Method, scalar version
// template<typename CallableF, typename CallableDF, typename InputType,
//          typename = std::enable_if_t<std::is_floating_point<InputType>::value>>
// InputType newton(CallableF&& f, CallableDF&& df, InputType x, double tol=1.0e-6, std::size_t maxIter=20){
//     static_assert(std::is_floating_point<decltype(f(x))>::value, "Function must return a scalar.");
//     static_assert(std::is_floating_point<decltype(df(x))>::value, "Derivative function must return a scalar.");
//     for (std::size_t iter = 0; iter < maxIter; iter++){
//         InputType dfx = df(x);
//         if (dfx == 0.0) throw NewtonSingularityError("ERROR: Singularity encountered. Newton's method cannot continue.");
//         InputType dx = -f(x)/dfx;
//         x += dx;
//         if (x == 0.0){
//             if (std::abs(dx) < tol) return x;
//         } else if (std::abs(dx/x) < tol) return x;
//     }
//     throw NewtonMaxIterationsError("ERROR: Maximum number of iterations exceeded.");
// }

// template<typename Callable, typename InputType,
//          typename = std::enable_if_t<std::is_floating_point<InputType>::value>>
// InputType newton(Callable&& f, InputType x, double tol=1.0e-6, std::size_t maxIter=20){
//     return newton(std::forward<Callable>(f),
//                   [&f](InputType x0){ return df(std::forward<Callable>(f), x0); },
//                   x, tol, maxIter);
// }

// Newton's method, vector version
template<typename CallableF, typename CallableJ, typename InputType,
         typename = std::enable_if_t<std::decay_t<InputType>::RowsAtCompileTime == 1
                                    || std::decay_t<InputType>::ColsAtCompileTime == 1>>
std::decay_t<InputType> newton(CallableF&& f, CallableJ&& J, InputType&& x, double tol=1.0e-6, std::size_t maxIter=20){
    bool dimensionsChecked = false;
    for (std::size_t iter = 0; iter < maxIter; iter++){
        auto fx = f(std::forward<InputType>(x));
        auto Jx = J(std::forward<InputType>(x));
        std::cout << "iter = " << iter << ": x = " << x.transpose() << ", f(x) = " << fx.transpose() << std::endl;
        std::cout << "Jx = " << std::endl << Jx << std::endl;
        if (!dimensionsChecked){
            // Check against dimension mismatch, for the first iteration only
            if (fx.size() < x.size() || fx.size() < Jx.rows())
                throw std::invalid_argument("ERROR: Output dimension must not be less than input dimension.");
            if (fx.size() != Jx.cols()) throw std::invalid_argument("ERROR: Function and derivative must have the same dimension.");
            dimensionsChecked = true;
        }
        Eigen::ColPivHouseholderQR<Eigen::Ref<std::decay_t<decltype(Jx)> > > qr(Jx);
        if (qr.rank() < fx.size()) throw NewtonSingularityError("ERROR: Singularity encountered. Newton's method cannot continue.");
        auto dx = qr.solve(fx);
        x -= dx;

        std::cout << "Error = " << dx.squaredNorm() / x.squaredNorm() << std::endl << std::endl;
        if (x.squaredNorm() == 0.0) {
            if (dx.squaredNorm() < tol*tol) return x;
        } else if (dx.squaredNorm() / x.squaredNorm() < tol*tol) return x;
    }
    throw NewtonMaxIterationsError("ERROR: Maximum number of iterations exceeded.");
}

template<typename Callable, typename InputType,
         typename = std::enable_if_t<std::decay_t<InputType>::RowsAtCompileTime == 1
                                    || std::decay_t<InputType>::ColsAtCompileTime == 1>>
std::decay_t<InputType> newton(Callable&& f, InputType&& x, double tol=1.0e-6, std::size_t maxIter=20){
    return newton(std::forward<Callable>(f),
                  [&f](InputType x0){ return Jacobian(std::forward<Callable>(f), x0); },
                  std::forward<InputType>(x), tol, maxIter);
}