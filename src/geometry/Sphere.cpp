#include "Sphere.h"
#include "Box.h"
#include "Direction.h"
#include <cmath>
#include <exception>

static const double PI = acos(-1.0);

Sphere::Sphere(const Point& pt, double R, double Sigma_t, double refrac):
    Shape(Sigma_t, refrac),
    _origin(pt),
    _radius(R)
{
    if(_radius <= 0) throw std::invalid_argument("ERROR: Radius must be positive.");
}

Sphere::Sphere(double x, double y, double z, double R, double Sigma_t, double refrac):
    Shape(Sigma_t, refrac),
    _origin(x,y,z),
    _radius(R)
{
    if(_radius <= 0) throw std::invalid_argument("ERROR: Radius must be positive.");
}

void Sphere::setRadius(double R){
    if (R <= 0.0) throw std::invalid_argument("ERROR: Radius must be positive.");
    _radius = R;
}

double Sphere::surfaceArea() const noexcept{ return 4.0*PI*_radius*_radius; }
double Sphere::volume() const noexcept{ return 4.0/3*PI*_radius*_radius*_radius; }

bool Sphere::surfaceContains(const Point& p) const noexcept{
    double dx = _origin.x() - p.x();
    double dy = _origin.y() - p.y();
    double dz = _origin.z() - p.z();
    return fabs(dx*dx + dy*dy + dz*dz - _radius*_radius) <= Shape::eps;
}

bool Sphere::encloses(const Point& p) const noexcept{
    double dx = _origin.x() - p.x();
    double dy = _origin.y() - p.y();
    double dz = _origin.z() - p.z();
    return dx*dx + dy*dy + dz*dz - _radius*_radius < -Shape::eps;
}

bool Sphere::encloses(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    if (const Sphere* shape = dynamic_cast<const Sphere*>(&other)){
        if (shape->_radius > _radius) return false;
        if (shape->_radius == _radius) return shape->_origin == _origin;
        double dx = _origin.x() - shape->_origin.x();
        double dy = _origin.y() - shape->_origin.y();
        double dz = _origin.z() - shape->_origin.z();
        double diffRadii = _radius - shape->_radius;
        return dx*dx + dy*dy + dz*dz < diffRadii*diffRadii;         
    } else{ // Shape is a Box
        auto squared = [](double x) -> double { return x*x; };
        double dx2 = std::max(squared(_origin.x() - other.xMin()), squared(_origin.x() - other.xMax()));
        double dy2 = std::max(squared(_origin.y() - other.yMin()), squared(_origin.y() - other.yMax()));
        double dz2 = std::max(squared(_origin.z() - other.zMin()), squared(_origin.z() - other.zMax()));
        return dx2 + dy2 + dz2 < _radius*_radius;
    }
}

bool Sphere::overlaps(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    if (const Sphere* shape = dynamic_cast<const Sphere*>(&other)){
        double dx = _origin.x() - shape->_origin.x();
        double dy = _origin.y() - shape->_origin.y();
        double dz = _origin.z() - shape->_origin.z();
        double sumRadii = _radius + shape->_radius;
        return dx*dx + dy*dy + dz*dz < sumRadii*sumRadii;
    } else{ // Shape is a Box
        double dx = _origin.x() - std::max(other.xMin(), std::min(_origin.x(), other.xMax()));
        double dy = _origin.y() - std::max(other.yMin(), std::min(_origin.y(), other.yMax()));
        double dz = _origin.z() - std::max(other.zMin(), std::min(_origin.z(), other.zMax()));
        return dx*dx + dy*dy + dz*dz < _radius*_radius;        
    }
}

double Sphere::distanceToSurface(const Point& p, const Direction& dir) const noexcept{
    double dx = p.x() - _origin.x();
    double dy = p.y() - _origin.y();
    double dz = p.z() - _origin.z();
    double C = dx*dx + dy*dy + dz*dz - _radius*_radius;
    double B = 2*(dx*dir.dx() + dy*dir.dy() + dz*dir.dz());

    if (fabs(C) < Shape::eps){ // on the surface
        if (B >= 0.0) return 0.0; // leaves the surface
        return -B; // enters the surface and leaves at another location
    }

    double discr = B*B - 4*C;
    if (discr < 0) return NAN; // particle is outside and will never enters the Sphere
    double sp, sm; // two roots of the solutions
    if (B < 0){
        sp = (-B + sqrt(discr))/2;
        sm = C/sp; 
    } else{
        sm = (-B - sqrt(discr))/2;
        sp = C/sm;
    }

    if (C > 0) return (sp > 0 ? std::min(sp,sm) : NAN); // particle is outside, the roots have the same sign
    return (sp > 0 ? sp : sm); // the roots have opposite sign, take the positive one - particle is inside
}

Direction Sphere::normal(const Point& pos) const{
    if (!surfaceContains(pos)) throw std::invalid_argument("ERROR: Point is not on the surface.");
    return Direction(_origin, pos);
}

std::ostream& Sphere::print(std::ostream& os) const noexcept{
    if (_origin.x() == 0.0) os << "x^2";
    else{
        os << "(x ";
        if (_origin.x() > 0) os << "- " << _origin.x() << ")^2";
        else os << "+ " << -_origin.x() << ")^2";
    }
    os << " + ";

    if (_origin.y() == 0.0) os << "y^2";
    else{
        os << "(y ";
        if (_origin.y() > 0) os << "- " << _origin.y() << ")^2";
        else os << "+ " << -_origin.y() << ")^2";
    }
    os << " + ";

    if (_origin.z() == 0.0) os << "z^2";
    else{
        os << "(z ";
        if (_origin.z() > 0) os << "- " << _origin.z() << ")^2";
        else os << "+ " << -_origin.z() << ")^2";
    }
    os << " = " << _radius*_radius;
    return os;
}