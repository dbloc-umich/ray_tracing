#include "Box.h"
#include "Direction.h"
#include <cmath>
#include <exception>

Box::Box(const Point& lower, const Point& upper, double Sigma_t, double refrac):
    Shape(Sigma_t, refrac),
    _lower(std::min(lower.x(), upper.x()), std::min(lower.y(), upper.y()), std::min(lower.z(), upper.z())),
    _upper(std::max(lower.x(), upper.x()), std::max(lower.y(), upper.y()), std::max(lower.z(), upper.z()))
{
    if(_lower.x() == _upper.x() || _lower.y() == _upper.y() || _lower.z() == _upper.z())
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
}

Box::Box(double x0, double y0, double z0, double x1, double y1, double z1, double Sigma_t, double refrac):
    Shape(Sigma_t, refrac),
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
    if (fabs(p.x()-xMin()) <= Shape::eps || fabs(p.x()-xMax()) <= Shape::eps)
        return p.y() > yMin() && p.y() < yMax() && p.z() > zMin() && p.z() < zMax();
    if (fabs(p.y()-yMin()) <= Shape::eps || fabs(p.y()-yMax()) <= Shape::eps)
        return p.z() > zMin() && p.z() < zMax() && p.x() > xMin() && p.x() < xMax();
    if (fabs(p.z()-zMin()) <= Shape::eps || fabs(p.z()-zMax()) <= Shape::eps)
        return p.x() > xMin() && p.x() < xMax() && p.y() > yMin() && p.y() < yMax();
    return false;
}

bool Box::encloses(const Point& p) const noexcept{
    // If this->surfaceContains(p) is true, then this->encloses(p) is false;
    bool x = xMin()-p.x() < -Shape::eps && xMax()-p.x() > Shape::eps;
    bool y = yMin()-p.y() < -Shape::eps && yMax()-p.y() > Shape::eps;
    bool z = zMin()-p.z() < -Shape::eps && zMax()-p.z() > Shape::eps;
    return x && y && z;
}

bool Box::encloses(const Shape& other) const noexcept{
    bool x = xMin() < other.xMin() && xMax() > other.xMax();
    bool y = yMin() < other.yMin() && yMax() > other.yMax();
    bool z = zMin() < other.zMin() && zMax() > other.zMax();
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
        if ((fabs(p.x()-xMin()) <= Shape::eps || fabs(p.x()-xMax()) <= Shape::eps) && dir.dx() == 0.0) return 0.0;
        if ((fabs(p.y()-yMin()) <= Shape::eps || fabs(p.y()-yMax()) <= Shape::eps) && dir.dy() == 0.0) return 0.0;
        if ((fabs(p.z()-zMin()) <= Shape::eps || fabs(p.z()-zMax()) <= Shape::eps) && dir.dz() == 0.0) return 0.0;
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
        if (!isOnSurface || fabs(p.x()-xMin()) > Shape::eps){ // Point not on the xMin surface
            d = (xMin()-p.x()) / dir.dx();
            if (isValidDistY(d) && isValidDistZ(d)){ s = (d < s && d > 0) ? d : s; }
        }
        if (!isOnSurface || fabs(p.x()-xMax()) > Shape::eps){ // Point not on the xMax surface
            d = (xMax()-p.x()) / dir.dx();
            if (isValidDistY(d) && isValidDistZ(d)){ s = (d < s && d > 0) ? d : s; }
        }
    }

    if (dir.dy() != 0.0){
        if (!isOnSurface || fabs(p.y()-yMin()) > Shape::eps){ // Point not on the yMin surface
            d = (yMin()-p.y()) / dir.dy();
            if (isValidDistZ(d) && isValidDistX(d)){ s = (d < s && d > 0) ? d : s; }
        }
        if (!isOnSurface || fabs(p.y()-yMax()) > Shape::eps){ // Point not on the yMax surface
            d = (yMax()-p.y()) / dir.dy();
            if (isValidDistZ(d) && isValidDistX(d)){ s = (d < s && d > 0) ? d : s; }
        }
    }

    if (dir.dz() != 0.0){
        if (!isOnSurface || fabs(p.z()-zMin()) > Shape::eps){ // Point not on the zMin surface
            d = (zMin()-p.z()) / dir.dz();
            if (isValidDistX(d) && isValidDistY(d)){ s = (d < s && d > 0) ? d : s; }
        }
        if (!isOnSurface || fabs(p.z()-zMax()) > Shape::eps){ // Point not on the zMax surface
            d = (zMax()-p.z()) / dir.dz();
            if (isValidDistX(d) && isValidDistY(d)){ s = (d < s && d > 0) ? d : s; }
        }
    }
    return s == std::numeric_limits<double>::max() ? (isOnSurface ? 0.0 : NAN) : s;
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