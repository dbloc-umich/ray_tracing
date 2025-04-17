#ifndef DIRECTION_H
#define DIRECTION_H
#include "Vector.h"

class Point;
class Direction: public Vector{
    public:
    explicit Direction(std::nullptr_t);
    Direction(double mu, double gamma);
    Direction(double dx, double dy, double dz);
    Direction(const Point& p1, const Point& p2);
    Direction(const Vector& v): Direction(v.dx(), v.dy(), v.dz()){}

    double norm() const noexcept override;
};

#endif // DIRECTION_H