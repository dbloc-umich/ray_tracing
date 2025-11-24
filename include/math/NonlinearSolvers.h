#include "Eigen/Dense"
#include "Derivative.h"
#include <cmath>
#include <type_traits>

// Failure cases
class NewtonSingularityError : public std::runtime_error {
public:
    explicit NewtonSingularityError(const std::string& msg) : std::runtime_error(msg) {}
};


// Thrown when Newton's method exceeds maximum allowed iterations
class NewtonMaxIterationsError : public std::runtime_error {
public:
    explicit NewtonMaxIterationsError(const std::string& msg) : std::runtime_error(msg) {}
};

// Newton's Method, scalar version
template<typename CallableF, typename CallableDF, typename InputType,
         typename = std::enable_if_t<std::is_floating_point<InputType>::value>>
InputType newton(CallableF&& f, CallableDF&& df, InputType x, double tol=1.0e-6, std::size_t maxIter=20){
    static_assert(std::is_floating_point<decltype(f(x))>::value, "Function must return a scalar.");
    static_assert(std::is_floating_point<decltype(df(x))>::value, "Derivative function must return a scalar.");
    for (std::size_t iter = 0; iter < maxIter; iter++){
        InputType dfx = df(x);
        if (dfx == 0.0) throw NewtonSingularityError("ERROR: Singularity encountered. Newton's method cannot continue.");
        InputType dx = -f(x)/dfx;
        x += dx;
        if (x == 0.0){
            if (std::abs(dx) < tol) return x;
        } else if (std::abs(dx/x) < tol) return x;
    }
    throw NewtonMaxIterationsError("ERROR: Maximum number of iterations exceeded.");
}

template<typename Callable, typename InputType,
         typename = std::enable_if_t<std::is_floating_point<InputType>::value>>
InputType newton(Callable&& f, InputType x, double tol=1.0e-6, std::size_t maxIter=20){
    return newton(std::forward<Callable>(f),
                  [&f](InputType x0){ return df(std::forward<Callable>(f), x0); },
                  x, tol, maxIter);
}

// Newton's method, vector version
template<typename CallableF, typename CallableJ, typename InputType,
         typename = std::enable_if_t<std::decay_t<InputType>::RowsAtCompileTime == 1
                                    || std::decay_t<InputType>::ColsAtCompileTime == 1>>
std::decay_t<InputType> newton(CallableF&& f, CallableJ&& J, InputType&& x, double tol=1.0e-6, std::size_t maxIter=20){
    bool dimensionsChecked = false;
    for (std::size_t iter = 0; iter < maxIter; iter++){
        auto fx = f(std::forward<InputType>(x));
        auto Jx = J(std::forward<InputType>(x));
        if (!dimensionsChecked){
            // Check against dimension mismatch, for the first iteration only
            if (fx.size() != x.size()) throw std::invalid_argument("ERROR: Input and output dimensions must match.");
            if (fx.size() != Jx.rows() || fx.size() != Jx.cols())
                throw std::invalid_argument("ERROR: Function and derivative must have the same dimension.");
            dimensionsChecked = true;
        }
        auto QrFactored = Jx.colPivHouseholderQr();
        if (QrFactored.rank() < Jx.cols()) throw NewtonSingularityError("ERROR: Singularity encountered. Newton's method cannot continue.");
        auto dx = -QrFactored.solve(fx);
        x += dx;
        if (x.squaredNorm() == 0.0) {
            if (dx.squaredNorm() < tol) return x;
        } else if (dx.squaredNorm() / x.squaredNorm() < tol) return x;
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