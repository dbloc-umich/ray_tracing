#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include "Eigen/Dense"

// First derivative
// Derivative of a single-variable function
template<typename Callable, typename Scalar,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto df(Callable&& f, Scalar x) -> decltype(f(x));

// Partial derivative of a multi-variable function, with respect to variable x[ind]
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto df(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x, Eigen::Index ind) -> decltype(f(x));

// First derivative operators
// Gradient
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto grad(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x);

// Divergent
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
Scalar div(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x);

// Curl
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto curl(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x);

// Jacobian
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto Jacobian(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x);

// Second derivative
// Derivative of a single-variable function
template<typename Callable, typename Scalar,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto df2(Callable&& f, Scalar x) -> decltype(f(x));

// Partial derivative of a multi-variable function, with respect to variable x[ind]
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto df2(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x, Eigen::Index ind1, Eigen::Index ind2) -> decltype(f(x));

// Hessian
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point_v<Scalar>>>
auto Hessian(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x);

#include "Derivative.tpp"
#endif