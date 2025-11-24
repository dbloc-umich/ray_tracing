#include "Sphere.h"
#include "Box.h"
#include "Constants.h"
#include "Ellipsoid.h"

#include <cmath>

Sphere::Sphere(const Eigen::Vector3d& pt, double R, std::shared_ptr<Material> mat):
    Shape(mat),
    _center(pt),
    _radius(R)
{
    if(_radius <= 0) throw std::invalid_argument("ERROR: Radius must be positive.");
}

Sphere::Sphere(double x, double y, double z, double R, std::shared_ptr<Material> mat):
    Shape(mat),
    _center(x,y,z),
    _radius(R)
{
    if(_radius <= 0) throw std::invalid_argument("ERROR: Radius must be positive.");
}

void Sphere::setRadius(double R){
    if (R <= 0.0) throw std::invalid_argument("ERROR: Radius must be positive.");
    _radius = R;
}

double Sphere::surfaceArea() const noexcept{ return 4.0*mconst::pi*_radius*_radius; }
double Sphere::volume() const noexcept{ return 4.0/3*mconst::pi*_radius*_radius*_radius; }

bool Sphere::surfaceContains(const Eigen::Vector3d& p) const noexcept{
    return std::abs((_center - p).squaredNorm() - _radius*_radius) <= Shape::eps;
}

bool Sphere::encloses(const Eigen::Vector3d& p) const noexcept{
    return (_center - p).squaredNorm() - _radius*_radius < -Shape::eps;
}

// bool Sphere::encloses(const Shape& other) const noexcept{
//     if (this == &other) return true;  // same object
//     if (const Sphere* shape = dynamic_cast<const Sphere*>(&other)){
//         if (shape->_radius > _radius) return false;
//         if (shape->_radius == _radius) return shape->_center == _center;
//         Eigen::Vector3d delr = _center - shape->_center;
//         double diffRadii = _radius - shape->_radius;
//         return delr.squaredNorm() < diffRadii*diffRadii;         
//     } else{ // Shape is a Box
//         auto squared = [](double x) -> double { return x*x; };
//         double dx2 = std::max(squared(_center.x() - other.xMin()), squared(_center.x() - other.xMax()));
//         double dy2 = std::max(squared(_center.y() - other.yMin()), squared(_center.y() - other.yMax()));
//         double dz2 = std::max(squared(_center.z() - other.zMin()), squared(_center.z() - other.zMax()));
//         return dx2 + dy2 + dz2 < _radius*_radius;
//     }
//     return true;
// }

bool Sphere::overlaps(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    if (const Sphere* shape = dynamic_cast<const Sphere*>(&other)){
        Eigen::Vector3d delr = _center - shape->_center;
        double R = _radius + shape->_radius;
        return delr.squaredNorm() < R*R;
    }
    if (dynamic_cast<const Box*>(&other)){
        Eigen::Vector3d delr;
        delr[0] = _center[0] - std::max(other.xMin(), std::min(_center[0], other.xMax()));
        delr[1] = _center[1] - std::max(other.yMin(), std::min(_center[1], other.yMax()));
        delr[2] = _center[2] - std::max(other.zMin(), std::min(_center[2], other.zMax()));
        return delr.squaredNorm() < _radius*_radius;
    }
    return other.encloses(*this); // See the implementation files for other Shapes
}

double Sphere::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{

    Eigen::Vector3d delr = p - _center;
    double C = delr.squaredNorm() - _radius*_radius;
    double B = delr.dot(dir.value());

    if (std::abs(C) < Shape::eps){ // on the surface
        if (B >= 0.0) return 0.0; // leaves the surface
        return -2*B; // enters the surface and leaves at another location
    }
    
    double discr = B*B - C;
    if (discr < 0) return NAN; // particle is outside and will never enters the Sphere
    double sp, sm; // two roots of the solutions
    if (B < 0){
        sp = -B + std::sqrt(discr);
        sm = C/sp; 
    } else{
        sm = -B - std::sqrt(discr);
        sp = C/sm;
    }

    if (C > 0) return (sp > 0 ? std::min(sp,sm) : NAN); // particle is outside, the roots have the same sign
    return (sp > 0 ? sp : sm); // the roots have opposite sign, take the positive one - particle is inside
}

UnitVector3d Sphere::normal(const Eigen::Vector3d& pos) const{
    if (!surfaceContains(pos)) throw std::invalid_argument("ERROR: Eigen::Vector3d is not on the surface.");
    return UnitVector3d(pos - _center);
}

std::ostream& Sphere::print(std::ostream& os) const noexcept{
    if (_center[0] == 0.0) os << "x^2";
    else{
        os << "(x ";
        if (_center.x() > 0) os << "- " << _center[0] << ")^2";
        else os << "+ " << -_center[0] << ")^2";
    }
    os << " + ";

    if (_center[1] == 0.0) os << "y^2";
    else{
        os << "(y ";
        if (_center[1] > 0) os << "- " << _center[1] << ")^2";
        else os << "+ " << -_center[1] << ")^2";
    }
    os << " + ";

    if (_center[2] == 0.0) os << "z^2";
    else{
        os << "(z ";
        if (_center.z() > 0) os << "- " << _center[2] << ")^2";
        else os << "+ " << -_center[2] << ")^2";
    }
    os << " = " << _radius*_radius;
    return os;
}