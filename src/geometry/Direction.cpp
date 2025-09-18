#include "Direction.h"
#include "Point.h"

#include <cmath>

Direction::Direction(std::nullptr_t):
    Vector(NAN, NAN, NAN)
{}

Direction::Direction(double mu, double gamma):
    Vector( fabs(mu) > 1.0 ? NAN : mu, std::sqrt(1-mu*mu)*cos(gamma), std::sqrt(1-mu*mu)*sin(gamma))
{}

Direction::Direction(double dx, double dy, double dz):
    Vector(dx, dy, dz)
{
    *this /= Vector::norm();
}

Direction::Direction(const Point& p1, const Point& p2):
    Vector(p1, p2)
{
    *this /= Vector::norm();
}

double Direction::norm() const noexcept{ return *this ? 1.0 : NAN; }