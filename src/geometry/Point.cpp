#include "Point.h"
#include "Vector.h"

#include <cmath>
#include <limits>

static constexpr double eps = 1e-9;
double Point::dist(const Point& other) const noexcept{
    return sqrt((other._x-_x)*(other._x-_x) + (other._y-_y)*(other._y-_y) + (other._z-_z)*(other._z-_z));
}

bool Point::operator==(const Point& other) const noexcept{
    bool x = fabs(_x - other._x) <= eps;
    bool y = fabs(_y - other._y) <= eps;
    bool z = fabs(_z - other._z) <= eps;
    return x && y && z;
}

Point& Point::advance(const Vector& v, const double s){
    _x += v.dx()*s;
    _y += v.dy()*s;
    _z += v.dz()*s;
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Point& p){
    os << "(" << p._x << ", " << p._y << ", " << p._z << ")";
    return os;
}