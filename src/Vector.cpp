#include "Point.h"
#include "Vector.h"

#include <cmath>
#include <limits>

static constexpr double eps = 1e-9;
Vector::Vector(double dx, double dy, double dz):
    _dx(dx),
    _dy(dy),
    _dz(dz)
{}

Vector::Vector(const Point& p1, const Point& p2):
    _dx(p2.x()-p1.x()),
    _dy(p2.y()-p1.y()),
    _dz(p2.z()-p1.z())
{}

Vector Vector::operator+() const noexcept{ return Vector(*this); }
Vector Vector::operator+(const Vector& other) const noexcept{ return Vector(_dx+other._dx, _dy+other._dy, _dz+other._dz); }
Vector& Vector::operator+=(const Vector& other) noexcept{
    _dx += other._dx;
    _dy += other._dy;
    _dz += other._dz;
    return *this;
}

Vector Vector::operator-() const noexcept{ return Vector(-_dx, -_dy, -_dz); }
Vector Vector::operator-(const Vector& other) const noexcept{ return Vector(_dx-other._dx, _dy-other._dy, _dz-other._dz); }
Vector& Vector::operator-=(const Vector& other) noexcept{
    _dx -= other._dx;
    _dy -= other._dy;
    _dz -= other._dz;
    return *this;
}

Vector Vector::operator*(double d) const noexcept{ return Vector(_dx*d, _dy*d, _dz*d); }
Vector& Vector::operator*=(double d) noexcept{
    _dx *= d;
    _dy *= d;
    _dz *= d;
    return *this;
}

Vector Vector::operator/(double d) const noexcept{ return Vector(_dx/d, _dy/d, _dz/d); }
Vector& Vector::operator/=(double d) noexcept{
    _dx /= d;
    _dy /= d;
    _dz /= d;
    return *this;
}

bool Vector::operator==(const Vector& other) const noexcept{
    bool x = fabs(_dx-other._dx) <= eps;
    bool y = fabs(_dy-other._dy) <= eps;
    bool z = fabs(_dz-other._dz) <= eps;
    return x && y && z;
}

Vector::operator bool() const noexcept{ return !(std::isnan(_dx) || std::isnan(_dy) || std::isnan(_dz)); }

double Vector::dot(const Vector& other) const noexcept{ return _dx*other._dx + _dy*other._dy + _dz*other._dz; }
double Vector::norm() const noexcept{ return sqrt(this->dot(*this)); }
Vector Vector::cross(const Vector& other) const noexcept{
    double dx = _dy*other._dz - _dz*other._dy;
    double dy = _dz*other._dx - _dx*other._dz;
    double dz = _dx*other._dy - _dy*other._dx;
    return Vector(dx, dy, dz);
}

bool Vector::isOrthogonal(const Vector& other) const noexcept{ return fabs(this->dot(other)) <= eps; }
bool Vector::isParallel(const Vector& other) const noexcept{
    Vector v = this->cross(other);
    return fabs(v._dx) <= eps && fabs(v._dy) <= eps && fabs(v._dz) <= eps; 
}

std::ostream& operator<<(std::ostream& os, const Vector& v){
    os << "[" << v._dx << ", " << v._dy << ", " << v._dz << "]";
    return os;
}