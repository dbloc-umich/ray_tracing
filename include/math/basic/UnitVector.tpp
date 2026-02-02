#include "UnitVector.h"

template<typename Scalar, int n>
UnitVector<Scalar, n>::UnitVector(const Eigen::Matrix<Scalar, n, 1>& v):
    _v(v)
{
    double norm = _v.squaredNorm();
    if (_v.squaredNorm() < eps) _v.fill(std::numeric_limits<Scalar>::quiet_NaN());
    else{
        norm = std::sqrt(norm);
        _v /= norm;
    }
}

template<typename Scalar, int n>
UnitVector<Scalar, n> UnitVector<Scalar, n>::operator-() const noexcept{
    UnitVector<Scalar, n> other = *this;
    other._v *= -1;
    return other;
}

template<typename Scalar, int n>
UnitVector<Scalar, n> UnitVector<Scalar, n>::reflect(const UnitVector<Scalar, n>& normal) const{
    auto res = _v - 2*(normal._v.dot(_v))*normal._v;
    return UnitVector<Scalar, n>(res);  
}

template<typename Scalar, int n>
UnitVector<Scalar, 3> UnitVector<Scalar, n>::refract(const UnitVector<Scalar, 3>& normal, double n1, double n2) const{
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
std::ostream& operator<<(std::ostream& os, const UnitVector<T, N>& v){
    os << "[" << v.value().transpose() << "]";
    return os;
}