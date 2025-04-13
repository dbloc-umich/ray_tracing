#include "Direction.h"
#include "Point.h"

#include <cmath>
#include <exception>

Direction::Direction(double mu, double gamma):
    Vector(mu, sqrt(1-mu*mu)*cos(gamma), sqrt(1-mu*mu)*sin(gamma))
{
    if(mu < -1.0 || mu > 1.0) throw std::domain_error("ERROR: The polar direction must be within the domain [-1,1].");
}

Direction::Direction(double dx, double dy, double dz):
    Vector(dx, dy, dz)
{
    double norm = Vector::norm();
    if (norm > 0){
        _dx /= norm;
        _dy /= norm;
        _dz /= norm;
    }
}

Direction::Direction(const Point& p1, const Point& p2):
    Vector(p1, p2)
{
    double norm = Vector::norm();
    if (norm > 0){
        _dx /= norm;
        _dy /= norm;
        _dz /= norm;
    }
}