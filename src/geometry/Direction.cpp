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

Direction reflected(const Direction& in, const Direction& normal){ return in - 2*(normal.dot(in))*normal; }
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2){
    if (in.isParallel(normal) || in.isOrthogonal(normal) || n1 == n2) return in;

    Vector c = in.cross(normal)*(n1/n2);
    double sine = c.norm();
    if (sine <= 1.0){
        Direction ortho = normal.cross(c); // component of the refracted vector that is orthogonal to normal
        int sgn = in.dot(normal) > 0.0 ? 1 : -1;
        return ortho*sine + sgn*normal*std::sqrt(1-sine*sine);
    }
    return Direction(nullptr);
}