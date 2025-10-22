#include "Box.h"
#include "Direction.h"
#include <cmath>
#include <limits>

Box::Box(const Point& lower, const Point& upper, const std::shared_ptr<Material>& mat):
    Shape(mat),
    _lower(std::min(lower.x(), upper.x()), std::min(lower.y(), upper.y()), std::min(lower.z(), upper.z())),
    _upper(std::max(lower.x(), upper.x()), std::max(lower.y(), upper.y()), std::max(lower.z(), upper.z()))
{
    if(_lower.x() == _upper.x() || _lower.y() == _upper.y() || _lower.z() == _upper.z())
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
}

Box::Box(double x0, double y0, double z0, double x1, double y1, double z1, const std::shared_ptr<Material>& mat):
    Shape(mat),
    _lower(std::min(x0,x1), std::min(y0,y1), std::min(z0,z1)),
    _upper(std::max(x0,x1), std::max(y0,y1), std::max(z0,z1))
{
    if(_lower.x() == _upper.x() || _lower.y() == _upper.y() || _lower.z() == _upper.z())
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
}

void Box::setLowerVertex(const Point& point){
    if (point.x() == _upper.x() || point.y() == _upper.y() || point.z() == _upper.z())
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
    _lower = point;
}

void Box::setUpperVertex(const Point& point){
    if (point.x() == _lower.x() || point.y() == _lower.y() || point.z() == _lower.z())
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate");
    _upper = point;
}

void Box::setVertices(const Point& lower, const Point& upper){
    if (lower.x() == upper.x() || lower.y() == upper.y() || lower.z() == upper.z())
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
    _lower = lower;
    _upper = upper;
}

bool Box::surfaceContains(const Point& p) const noexcept{
    bool inX = p.x()-xMin() >= -Shape::eps && p.x()-xMax() <= Shape::eps;
    bool inY = p.y()-yMin() >= -Shape::eps && p.y()-yMax() <= Shape::eps;
    bool inZ = p.z()-zMin() >= -Shape::eps && p.z()-zMax() <= Shape::eps;
    bool onX = fabs(p.x()-xMin()) <= Shape::eps || fabs(p.x()-xMax()) <= Shape::eps;
    bool onY = fabs(p.y()-yMin()) <= Shape::eps || fabs(p.y()-yMax()) <= Shape::eps;
    bool onZ = fabs(p.z()-zMin()) <= Shape::eps || fabs(p.z()-zMax()) <= Shape::eps;
    return (inX && inY && inZ) && (onX || onY || onZ);
}

bool Box::encloses(const Point& p) const noexcept{
    // If this->surfaceContains(p) is true, then this->encloses(p) is false;
    bool x = xMin()-p.x() < -Shape::eps && xMax()-p.x() > Shape::eps;
    bool y = yMin()-p.y() < -Shape::eps && yMax()-p.y() > Shape::eps;
    bool z = zMin()-p.z() < -Shape::eps && zMax()-p.z() > Shape::eps;
    return x && y && z;
}

bool Box::encloses(const Shape& other) const noexcept{
    bool x = xMin() <= other.xMin() && xMax() >= other.xMax();
    bool y = yMin() <= other.yMin() && yMax() >= other.yMax();
    bool z = zMin() <= other.zMin() && zMax() >= other.zMax();
    return x && y && z;
}

bool Box::overlaps(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    if (dynamic_cast<const Box*>(&other)){
        bool x = xMin() < other.xMax() && xMax() > other.xMin();
        bool y = yMin() < other.yMax() && yMax() > other.yMin();
        bool z = zMin() < other.zMax() && zMax() > other.zMin();
        return x && y && z;
    } else return other.overlaps(*this); // Shape is a Sphere
}

double Box::distanceToSurface(const Point& p, const Direction& dir) const noexcept{
    bool isOnSurface = surfaceContains(p);
    if (isOnSurface){ // point moving on a surface
        if (fabs(p.x()-xMin()) <= Shape::eps && dir.dx() <= 0.0) return 0.0;
        if (fabs(p.x()-xMax()) <= Shape::eps && dir.dx() >= 0.0) return 0.0;
        if (fabs(p.y()-yMin()) <= Shape::eps && dir.dy() <= 0.0) return 0.0;
        if (fabs(p.y()-yMax()) <= Shape::eps && dir.dy() >= 0.0) return 0.0;
        if (fabs(p.y()-zMin()) <= Shape::eps && dir.dz() <= 0.0) return 0.0;
        if (fabs(p.y()-zMax()) <= Shape::eps && dir.dz() >= 0.0) return 0.0;
    }
    
    double d;
    double s = std::numeric_limits<double>::max();
    auto distX = [&](double l) { return p.x() + dir.dx()*l; };
    auto distY = [&](double l) { return p.y() + dir.dy()*l; };
    auto distZ = [&](double l) { return p.z() + dir.dz()*l; };
    auto isValidDistX = [&](double l) { return distX(l) >= xMin() && distX(l) <= xMax(); };
    auto isValidDistY = [&](double l) { return distY(l) >= yMin() && distY(l) <= yMax(); };
    auto isValidDistZ = [&](double l) { return distZ(l) >= zMin() && distZ(l) <= zMax(); };

    if (dir.dx() != 0.0){
        double dest = NAN;
        if (p.x()-xMin() < -Shape::eps && dir.dx() > 0.0) dest = xMin(); // less than xMin
        else if (fabs(p.x()-xMin()) <= Shape::eps && dir.dx() > 0.0) dest = xMax(); // on xMin
        else if (p.x()-xMin() > Shape::eps && p.x()-xMax() < -Shape::eps) dest = dir.dx() > 0.0 ? xMax() : xMin(); // in between xMin and xMax
        else if (fabs(p.x()-xMax()) <= Shape::eps && dir.dx() < 0.0) dest = xMin(); // on xMax
        else if (p.x()-xMax() > Shape::eps && dir.dx() < 0.0) dest = xMax(); // greater than xMax

        if (std::isnan(dest)) return NAN;
        else{
            d = (dest-p.x())/dir.dx();
            if (isValidDistY(d) && isValidDistZ(d)){ s = (d < s) ? d : s; }
        }
    }

    if (dir.dy() != 0.0){
        double dest = NAN;
        if (p.y()-yMin() < -Shape::eps && dir.dy() > 0.0) dest = yMin(); // less than yMin
        else if (fabs(p.y()-yMin()) <= Shape::eps && dir.dy() > 0.0) dest = yMax(); // on yMin
        else if (p.y()-yMin() > Shape::eps && p.y()-xMax() < -Shape::eps) dest = dir.dy() > 0.0 ? yMax() : yMin(); // in between yMin and yMax
        else if (fabs(p.y()-yMax()) <= Shape::eps && dir.dy() < 0.0) dest = yMin(); // on yMax
        else if (p.y()-yMax() > Shape::eps && dir.dy() < 0.0) dest = yMax(); // greater than yMax

        if (std::isnan(dest)) return NAN;
        else{
            d = (dest-p.y())/dir.dy();
            if (isValidDistZ(d) && isValidDistX(d)){ s = (d < s) ? d : s; }
        }
    }

    if (dir.dz() != 0.0){
        double dest = NAN;
        if (p.z()-zMin() < -Shape::eps && dir.dz() > 0.0) dest = zMin(); // less than zMin
        else if (fabs(p.z()-zMin()) <= Shape::eps && dir.dz() > 0.0) dest = zMax(); // on zMin
        else if (p.z()-zMin() > Shape::eps && p.z()-zMax() < -Shape::eps) dest = dir.dz() > 0.0 ? zMax() : zMin(); // in between zMin and zMax
        else if (fabs(p.z()-zMax()) <= Shape::eps && dir.dz() < 0.0) dest = zMin(); // on zMax
        else if (p.z()-zMax() > Shape::eps && dir.dz() < 0.0) dest = zMax(); // greater than zMax

        if (std::isnan(dest)) return NAN;
        else{
            d = (dest-p.z())/dir.dz();
            if (isValidDistX(d) && isValidDistY(d)){ s = (d < s) ? d : s; }
        }
    }
    return s == std::numeric_limits<double>::max() ? NAN : s;
}

Point Box::centroid() const noexcept{
    double x = (xMin() + xMax())/2;
    double y = (yMin() + yMax())/2;
    double z = (zMin() + zMax())/2;
    return Point(x, y, z);
}

Direction Box::normal(const Point& pos) const{
    if (!surfaceContains(pos)) throw std::invalid_argument("ERROR: Point is not on the surface.");
    if (fabs(pos.x() - xMin()) < Shape::eps) return Direction(-1, 0, 0);
    if (fabs(pos.x() - xMax()) < Shape::eps) return Direction(1, 0, 0);
    if (fabs(pos.y() - yMin()) < Shape::eps) return Direction(0, -1, 0);
    if (fabs(pos.y() - yMax()) < Shape::eps) return Direction(0, 1, 0);
    if (fabs(pos.z() - zMin()) < Shape::eps) return Direction(0, 0, -1);
    return Direction(0, 0, 1);
}

std::ostream& Box::print(std::ostream& os) const noexcept{
    os << "[" << xMin() << ", " << xMax() << "] * ";
    os << "[" << yMin() << ", " << yMax() << "] * ";
    os << "[" << zMin() << ", " << zMax() << "]";
    return os;
}