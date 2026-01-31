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
    UnitVector(const Eigen::Matrix<Scalar, n, 1>& v):
        _v(v)
    {
        double norm = _v.squaredNorm();
        if (_v.squaredNorm() < eps) _v.fill(std::numeric_limits<Scalar>::quiet_NaN());
        else{
            norm = std::sqrt(norm);
            _v /= norm;
        }
    }

    template<typename First, typename... Rest,
            typename = std::enable_if_t<
                std::is_convertible_v<Scalar, std::decay_t<First>> &&
                (std::is_convertible_v<Scalar, std::decay_t<Rest>> && ...)
                >
            >
    explicit UnitVector(First&& first, Rest&&... rest)
        : UnitVector(Eigen::Matrix<Scalar, n, 1>{static_cast<Scalar>(std::forward<First>(first)),
                                                 static_cast<Scalar>(std::forward<Rest>(rest))...})
    {
        static_assert(sizeof...(Rest) + 1 == n, "Number of arguments must match the dimension of the vector.");
    }

    Eigen::Matrix<Scalar, n, 1>& value() noexcept{ return _v; }
    const Eigen::Matrix<Scalar, n, 1>& value() const noexcept{ return _v; }
    double operator[](Eigen::Index i) const noexcept{ return _v[i]; }
    UnitVector<Scalar, n> operator-() const noexcept{
        UnitVector<Scalar, n> other = *this;
        other._v *= -1;
        return other;
    }
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

    UnitVector<Scalar, n> reflect(const UnitVector<Scalar, n>& normal) const{
        auto res = _v - 2*(normal._v.dot(_v))*normal._v;
        return UnitVector<Scalar, n>(res);
    }

    UnitVector<Scalar, 3> refract(const UnitVector<Scalar, 3>& normal, double n1, double n2) const{
        if (_v.size() != 2 && _v.size() != 3){
            throw std::invalid_argument("ERROR: Refraction is only physically meaningful with 2- and 3-D vectors.");
        }
        if (_v.size() == 2){
            UnitVector<Scalar, 3> v2(_v[0], _v[1], 0); // recast to 3D vector
            return v2.refract(normal, n1, n2);
        }

        if (n1 == n2 || this->isOrthogonal(normal) || this->isParallel(normal)) return *this;
        Eigen::Matrix<Scalar, 3, 1> c = _v.cross(normal._v)*(n1/n2);
        double cnorm = c.norm();
        if (cnorm <= 1.0){
            Eigen::Matrix<Scalar, 3, 1> ortho = normal._v.cross(c/cnorm); // component of the refracted vector that is orthogonal to normal
            int sgn = _v.dot(normal._v) > 0.0 ? 1 : -1;
            ortho *= cnorm;
            ortho += sgn*std::sqrt(1-cnorm*cnorm)*normal._v;
            return UnitVector<Scalar, 3>(ortho);
        }
        return UnitVector<Scalar, 3>(nullptr);
    }

    template<typename T, int N>
    friend std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v);

    protected:
    static constexpr double eps = 1e-9;
    Eigen::Matrix<Scalar, n, 1> _v;
};

template<typename T, int N>
std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v){
    os << "[" << v.value().transpose() << "]";
    return os;
}

using UnitVector3d = UnitVector<double, 3>;

#endif // VECTOR_H