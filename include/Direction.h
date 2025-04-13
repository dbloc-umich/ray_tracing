#ifndef DIRECTION_H
#define DIRECTION_H
#include "Vector.h"

class Point;
class Direction: public Vector{
    public:
    explicit Direction(void*): Vector(0.0, 0.0, 0.0){}
    Direction(double mu, double gamma);
    Direction(double dx, double dy, double dz);
    Direction(const Point& p1, const Point& p2);
    Direction(const Vector& v): Direction(v.dx(), v.dy(), v.dz()){}

    void setX(double) noexcept = delete;
    void setY(double) noexcept = delete;
    void setZ(double) noexcept = delete;

    explicit operator bool() const noexcept{ return _dx != 0.0 || _dy != 0.0 || _dz != 0.0; }
    double norm() const noexcept override{ return *this ? 1.0 : 0.0; }
};

#endif // DIRECTION_H