#include "UnitVector3d.h"

#include <cmath>

static constexpr double eps = 1e-9;
UnitVector3d::UnitVector3d(std::nullptr_t):
    _v{std::nan(""), std::nan(""), std::nan("")}
{}

UnitVector3d::UnitVector3d(const Eigen::Vector3d& v):
    _v(v)
{
    double norm = _v.squaredNorm();
    if (_v.squaredNorm() < eps){
        _v[0] = std::nan("");
        _v[1] = std::nan("");
        _v[2] = std::nan("");
    } else{
        norm = std::sqrt(norm);
        _v /= norm;
    }
}

UnitVector3d::operator bool() const noexcept{
    return !(std::isnan(_v[0]) || std::isnan(_v[1]) || std::isnan(_v[2]));
}

bool UnitVector3d::isOrthogonal(const UnitVector3d& other) const noexcept{ return fabs(_v.dot(other._v)) < eps; }
bool UnitVector3d::isParallel(const UnitVector3d& other) const noexcept{ return _v.cross(other._v).squaredNorm() < eps; }

UnitVector3d UnitVector3d::reflect(const UnitVector3d& normal){
    Eigen::Vector3d res = _v - 2*(normal._v.dot(_v))*normal._v;
    return UnitVector3d(res);
}

UnitVector3d UnitVector3d::refract(const UnitVector3d& normal, double n1, double n2){
    if (n1 == n2 || this->isOrthogonal(normal) || this->isParallel(normal)) return *this;

    Eigen::Vector3d c = _v.cross(normal._v)*(n1/n2);
    double cnorm = c.norm();
    if (cnorm <= 1.0){
        Eigen::Vector3d ortho = normal._v.cross(c/cnorm); // component of the refracted vector that is orthogonal to normal
        int sgn = _v.dot(normal._v) > 0.0 ? 1 : -1;
        ortho *= cnorm;
        ortho += sgn*std::sqrt(1-cnorm*cnorm)*normal._v;
        return UnitVector3d(ortho);
    }
    return UnitVector3d(nullptr);
}

std::ostream& operator<<(std::ostream& os, const UnitVector3d& v){
    os << "[" << v._v[0] << ", " << v._v[1] << ", " << v._v[2] << "]";
    return os;
}