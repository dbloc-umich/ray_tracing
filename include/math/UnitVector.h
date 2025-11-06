#ifndef UNIT_VECTOR_H
#define UNIT_VECTOR_H

#include "Eigen/Dense"
#include <cmath>
#include <iostream>

#if __cplusplus >= 201703L
#include <type_traits>
#else
#include "all_convertible_to.h"
#endif

template<typename type, int n>
class UnitVector{
    public:
    explicit UnitVector(std::nullptr_t){ _v.fill(std::nan("")); }
    UnitVector(const Eigen::Matrix<type, n, 1>& v):
        _v(v)
    {
        double norm = _v.squaredNorm();
        if (_v.squaredNorm() < eps) _v.fill(std::nan(""));
        else{
            norm = std::sqrt(norm);
            _v /= norm;
        }
    }

#if __cplusplus >= 201703L
    template<typename First, typename... Rest,
            typename = std::enable_if_t<
                std::is_convertible<std::decay_t<First>, type>::value &&
                (std::is_convertible<std::decay_t<Rest>, type>::value && ...)
                >
            >
#else
    template<typename First, typename... Rest,
            typename = std::enable_if_t<
                std::is_convertible<std::decay_t<First>, type>::value &&
                all_convertible_to<type, std::decay_t<Rest>...>::value
                >
            >    
#endif
    explicit UnitVector(First&& first, Rest&&... rest)
        : UnitVector(Eigen::Matrix<type, n, 1>{static_cast<type>(std::forward<First>(first)),
                                               static_cast<type>(std::forward<Rest>(rest))...})
    {
        static_assert(sizeof...(Rest) + 1 == n, "Number of arguments must match the dimension of the vector.");
    }

    Eigen::Matrix<type, n, 1>& value() noexcept{ return _v; }
    const Eigen::Matrix<type, n, 1>& value() const noexcept{ return _v; }
    double operator[](std::size_t i) const noexcept{ return _v[i]; }
    explicit operator Eigen::Matrix<type, n, 1>() const noexcept{ return _v; }
    explicit operator bool() const noexcept{ return !_v.array().isNaN().all(); }

    bool isOrthogonal(const UnitVector<type, n>& other) const{ return std::fabs(_v.dot(other._v)) < eps; }
    bool isParallel(const UnitVector<type, n>& other) const{ return _v.cross(other._v).squaredNorm() < eps; }

    UnitVector<type, n> reflect(const UnitVector<type, n>& normal) const{
        auto res = _v - 2*(normal._v.dot(_v))*normal._v;
        return UnitVector<type, n>(res);
    }

    UnitVector<type, 3> refract(const UnitVector<type, 3>& normal, double n1, double n2) const{
        if (_v.size() != 2 && _v.size() != 3){
            throw std::invalid_argument("ERROR: Refraction is only physically meaningful with 2- and 3-D vectors.");
        }
        if (_v.size() == 2){
            UnitVector<type, 3> v2(_v[0], _v[1], 0); // recast to 3D vector
            return v2.refract(normal, n1, n2);
        }

        if (n1 == n2 || this->isOrthogonal(normal) || this->isParallel(normal)) return *this;
        Eigen::Matrix<type, 3, 1> c = _v.cross(normal._v)*(n1/n2);
        double cnorm = c.norm();
        if (cnorm <= 1.0){
            Eigen::Matrix<type, 3, 1> ortho = normal._v.cross(c/cnorm); // component of the refracted vector that is orthogonal to normal
            int sgn = _v.dot(normal._v) > 0.0 ? 1 : -1;
            ortho *= cnorm;
            ortho += sgn*std::sqrt(1-cnorm*cnorm)*normal._v;
            return UnitVector<type, 3>(ortho);
        }
        return UnitVector<type, 3>(nullptr);
    }

    template<typename T, int N>
    friend std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v);

    protected:
    static constexpr double eps = 1e-9;
    Eigen::Matrix<type, n, 1> _v;
};

template<typename T, int N>
std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v){
    os << "[" << v.value().transpose() << "]";
    return os;
}

using UnitVector3d = UnitVector<double, 3>;

#endif // VECTOR_H