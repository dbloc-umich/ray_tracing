#ifndef UNIT_VECTOR_H
#define UNIT_VECTOR_H

#include "Eigen/Dense"
#include <cmath>
#include <iostream>
#include <type_traits>

template<typename Scalar, int n>
class UnitVector{
    public:
    explicit UnitVector(std::nullptr_t){ _v.fill(std::numeric_limits<Scalar>::quiet_NaN()); }
    UnitVector(const Eigen::Matrix<Scalar, n, 1>& v);

    template<typename First, typename... Rest,
            typename = std::enable_if_t<std::is_convertible_v<Scalar, std::decay_t<First>> &&
                                        (std::is_convertible_v<Scalar, std::decay_t<Rest>> && ...)>>
    explicit UnitVector(First&& first, Rest&&... rest):
        UnitVector(Eigen::Matrix<Scalar, n, 1>{static_cast<Scalar>(std::forward<First>(first)),
                                            static_cast<Scalar>(std::forward<Rest>(rest))...})
    {
        static_assert(sizeof...(Rest) + 1 == n, "Number of arguments must match the dimension of the vector.");
    }
    
    Eigen::Matrix<Scalar, n, 1>& value() noexcept{ return _v; }
    const Eigen::Matrix<Scalar, n, 1>& value() const noexcept{ return _v; }
    double operator[](Eigen::Index i) const noexcept{ return _v[i]; }
    UnitVector<Scalar, n> operator-() const noexcept;
    explicit operator Eigen::Matrix<Scalar, n, 1>() const noexcept{ return _v; }
    explicit operator bool() const noexcept{ return !_v.array().isNaN().all(); }

    bool operator==(const UnitVector<Scalar, n>& other) const noexcept{ return _v == other._v; }
    bool operator==(const Eigen::Matrix<Scalar, n, 1>& other) const noexcept{ return _v == other; }
    bool operator!=(const UnitVector<Scalar, n>& other) const noexcept{ return _v != other._v; }
    bool operator!=(const Eigen::Matrix<Scalar, n, 1>& other) const noexcept{ return _v != other; }
    bool isApprox(const UnitVector<Scalar, n>& other) const noexcept{ return _v.isApprox(other._v); }
    bool isApprox(const Eigen::Matrix<Scalar, n, 1>& other) const noexcept{ return _v.isApprox(other); }

    bool isOrthogonal(const UnitVector<Scalar, n>& other) const{ return std::abs(_v.dot(other._v)) < eps; }
    bool isParallel(const UnitVector<Scalar, n>& other) const{ return _v.cross(other._v).squaredNorm() < eps; }

    UnitVector<Scalar, n> reflect(const UnitVector<Scalar, n>& normal) const;
    UnitVector<Scalar, 3> refract(const UnitVector<Scalar, 3>& normal, double n1, double n2) const;

    template<typename T, int N>
    friend std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v);

    protected:
    static constexpr double eps = 1e-9;
    Eigen::Matrix<Scalar, n, 1> _v;
};

template<typename T, int N>
std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v);

using UnitVector3d = UnitVector<double, 3>;

#include "UnitVector.tpp"
#endif // UNIT_VECTOR_H