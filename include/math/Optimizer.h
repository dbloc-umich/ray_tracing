#include "Eigen/Cholesky"
#include "Eigen/Dense"
#include "Derivative.h"
#include <cmath>
#include <type_traits>

// Projected Newton's method to minimize an objective function, subject to a box constraint
template<typename CallableF, typename CallableG, typename CallableH, typename InputType,
         typename = std::enable_if_t<std::decay_t<InputType>::RowsAtCompileTime == 1
                                    || std::decay_t<InputType>::ColsAtCompileTime == 1>>
std::decay_t<InputType> projectedNewton(CallableF&& f, InputType& x, const InputType& L, const InputType& U,
                                        CallableG&& df, CallableH&& Hf, double tol=1.0e-6, std::size_t maxIter=30){
    static_assert(std::is_floating_point<decltype(f(x))>::value, "objective function must return a floating point value.");
    assert(df(x).size() == Hf(x).rows() && df(x).size() == Hf(x).cols()); //dimensions of the gradient and Hessian must match

    auto fx0 = f(x);
    constexpr double c = 1e-4;
    for (std::size_t iter = 0; iter < maxIter; iter++){
        auto grad = df(x);
        auto hess = Hf(x);

        // Check if the Hessian is positive definite
        Eigen::LLT<Eigen::Ref<std::decay_t<decltype(hess)>>> llt(hess);
        double lambda = 1.0e-8; // Levenberg–Marquardt damping constant
        while (llt.info() != Eigen::Success){
            hess = llt.matrixL();
            hess = hess * hess.transpose();
            for (Eigen::Index i = 0; i < grad.size(); i++) hess(i,i) += lambda;
            lambda *= 10;
            llt.compute(hess);
        }

        // Clamp the step size and advance
        std::decay_t<decltype(grad)> dx = -llt.solve(grad);
        double alpha = 1.0; // maximum step size before hitting a bound
        for (Eigen::Index i = 0; i < grad.size(); i++){
            if (dx[i] == 0.0) continue;
            alpha = dx[i] > 0 ? std::min(alpha, (U[i]-x[i])/dx[i]) : std::min(alpha, (L[i]-x[i])/dx[i]);
        }
        auto fx1 = f(x+alpha*dx);
        while (fx1 > fx0 + c*alpha*grad.dot(dx)){ // Armijo–Goldstein condition
            alpha *= 0.5;
            fx1 = f(x+alpha*dx);
        }
        x += alpha*dx;
 
        // Convergence check
        for (Eigen::Index i = 0; i < grad.size(); i++) dx[i] /= (x[i] == 0.0 ? 1 : x[i]);
        if (dx.squaredNorm() <= tol*tol || std::abs(1.0 - fx1/fx0) <= tol) return x;
        fx0 = fx1;
    }
    throw std::runtime_error("ERROR: Maximum number of iterations exceeded.");
}