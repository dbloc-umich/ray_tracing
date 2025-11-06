#include "Sphere.h"
#include "Box.h"
#include "Constants.h"

#include <cmath>

Sphere::Sphere(const Eigen::Vector3d& pt, double R, std::shared_ptr<Material> mat):
    Shape(mat),
    _origin(pt),
    _radius(R)
{
    if(_radius <= 0) throw std::invalid_argument("ERROR: Radius must be positive.");
}

Sphere::Sphere(double x, double y, double z, double R, std::shared_ptr<Material> mat):
    Shape(mat),
    _origin(x,y,z),
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
    return fabs((_origin - p).squaredNorm() - _radius*_radius) <= Shape::eps;
}

bool Sphere::encloses(const Eigen::Vector3d& p) const noexcept{
    return (_origin - p).squaredNorm() - _radius*_radius < -Shape::eps;
}

bool Sphere::encloses(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    if (const Sphere* shape = dynamic_cast<const Sphere*>(&other)){
        if (shape->_radius > _radius) return false;
        if (shape->_radius == _radius) return shape->_origin == _origin;
        Eigen::Vector3d delr = _origin - shape->_origin;
        double diffRadii = _radius - shape->_radius;
        return delr.squaredNorm() < diffRadii*diffRadii;         
    } else{ // Shape is a Box
        auto squared = [](double x) -> double { return x*x; };
        double dx2 = std::max(squared(_origin.x() - other.xMin()), squared(_origin.x() - other.xMax()));
        double dy2 = std::max(squared(_origin.y() - other.yMin()), squared(_origin.y() - other.yMax()));
        double dz2 = std::max(squared(_origin.z() - other.zMin()), squared(_origin.z() - other.zMax()));
        return dx2 + dy2 + dz2 < _radius*_radius;
    }
    return true;
}

bool Sphere::overlaps(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    
    Eigen::Vector3d delr;
    double R = _radius;
    if (const Sphere* shape = dynamic_cast<const Sphere*>(&other)){
        delr = _origin - shape->_origin;
        R += shape->_radius;
    } else{ // Shape is a Box
        delr[0] = _origin[0] - std::max(other.xMin(), std::min(_origin[0], other.xMax()));
        delr[1] = _origin[1] - std::max(other.yMin(), std::min(_origin[1], other.yMax()));
        delr[2] = _origin[2] - std::max(other.zMin(), std::min(_origin[2], other.zMax()));
    }
    return delr.squaredNorm() < R*R;
}

double Sphere::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{

    Eigen::Vector3d delr = p - _origin;
    double C = delr.squaredNorm() - _radius*_radius;
    double B = 2*(delr.dot(dir.value()));

    if (fabs(C) < Shape::eps){ // on the surface
        if (B >= 0.0) return 0.0; // leaves the surface
        return -B; // enters the surface and leaves at another location
    }

    double discr = B*B - 4*C;
    if (discr < 0) return NAN; // particle is outside and will never enters the Sphere
    double sp, sm; // two roots of the solutions
    if (B < 0){
        sp = (-B + std::sqrt(discr))/2;
        sm = C/sp; 
    } else{
        sm = (-B - std::sqrt(discr))/2;
        sp = C/sm;
    }

    if (C > 0) return (sp > 0 ? std::min(sp,sm) : NAN); // particle is outside, the roots have the same sign
    return (sp > 0 ? sp : sm); // the roots have opposite sign, take the positive one - particle is inside
}

UnitVector3d Sphere::normal(const Eigen::Vector3d& pos) const{
    if (!surfaceContains(pos)) throw std::invalid_argument("ERROR: Eigen::Vector3d is not on the surface.");
    return UnitVector3d(pos - _origin);
}

std::ostream& Sphere::print(std::ostream& os) const noexcept{
    if (_origin[0] == 0.0) os << "x^2";
    else{
        os << "(x ";
        if (_origin.x() > 0) os << "- " << _origin[0] << ")^2";
        else os << "+ " << -_origin[0] << ")^2";
    }
    os << " + ";

    if (_origin[1] == 0.0) os << "y^2";
    else{
        os << "(y ";
        if (_origin[1] > 0) os << "- " << _origin[1] << ")^2";
        else os << "+ " << -_origin[1] << ")^2";
    }
    os << " + ";

    if (_origin[2] == 0.0) os << "z^2";
    else{
        os << "(z ";
        if (_origin.z() > 0) os << "- " << _origin[2] << ")^2";
        else os << "+ " << -_origin[2] << ")^2";
    }
    os << " = " << _radius*_radius;
    return os;
}