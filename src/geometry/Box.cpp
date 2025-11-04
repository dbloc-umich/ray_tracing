#include "Box.h"
#include <cmath>
#include <limits>

Box::Box(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper, std::shared_ptr<Material> mat):
    Shape(mat),
    _lower(std::min(lower[0], upper[0]), std::min(lower[1], upper[1]), std::min(lower[2], upper[2])),
    _upper(std::max(lower[0], upper[0]), std::max(lower[1], upper[1]), std::max(lower[2], upper[2]))
{
    if(_lower[0] == _upper[0] || _lower[1] == _upper[1] || _lower[2] == _upper[2])
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
}

Box::Box(double x0, double y0, double z0, double x1, double y1, double z1, std::shared_ptr<Material> mat):
    Shape(mat),
    _lower(std::min(x0,x1), std::min(y0,y1), std::min(z0,z1)),
    _upper(std::max(x0,x1), std::max(y0,y1), std::max(z0,z1))
{
    if(_lower[0] == _upper[0] || _lower[1] == _upper[1] || _lower[2] == _upper[2])
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
}

void Box::setLowerVertex(const Eigen::Vector3d& point){
    if (point[0] == _upper[0] || point[1] == _upper[1] || point[2] == _upper[2])
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
    _lower = point;
}

void Box::setUpperVertex(const Eigen::Vector3d& point){
    if (point[0] == _lower[0] || point[1] == _lower[1] || point[2] == _lower[2])
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate");
    _upper = point;
}

void Box::setVertices(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper){
    if (lower[0] == upper[0] || lower[1] == upper[1] || lower[2] == upper[2])
        throw std::invalid_argument("ERROR: The vertices coincide in at least one coordinate.");
    _lower = lower;
    _upper = upper;
}

bool Box::surfaceContains(const Eigen::Vector3d& p) const noexcept{
    bool inX = p[0]-xMin() >= -Shape::eps && p[0]-xMax() <= Shape::eps;
    bool inY = p[1]-yMin() >= -Shape::eps && p[1]-yMax() <= Shape::eps;
    bool inZ = p[2]-zMin() >= -Shape::eps && p[2]-zMax() <= Shape::eps;
    bool onX = fabs(p[0]-xMin()) <= Shape::eps || fabs(p[0]-xMax()) <= Shape::eps;
    bool onY = fabs(p[1]-yMin()) <= Shape::eps || fabs(p[1]-yMax()) <= Shape::eps;
    bool onZ = fabs(p[2]-zMin()) <= Shape::eps || fabs(p[2]-zMax()) <= Shape::eps;
    return (inX && inY && inZ) && (onX || onY || onZ);
}

bool Box::encloses(const Eigen::Vector3d& p) const noexcept{
    // If this->surfaceContains(p) is true, then this->encloses(p) is false;
    bool x = xMin()-p[0] < -Shape::eps && xMax()-p[0] > Shape::eps;
    bool y = yMin()-p[1] < -Shape::eps && yMax()-p[1] > Shape::eps;
    bool z = zMin()-p[2] < -Shape::eps && zMax()-p[2] > Shape::eps;
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

double Box::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{
    bool isOnSurface = surfaceContains(p);
    if (isOnSurface){ // point moving on a surface
        if (fabs(p[0]-xMin()) <= Shape::eps && dir[0] <= 0.0) return 0.0;
        if (fabs(p[0]-xMax()) <= Shape::eps && dir[0] >= 0.0) return 0.0;
        if (fabs(p[1]-yMin()) <= Shape::eps && dir[1] <= 0.0) return 0.0;
        if (fabs(p[1]-yMax()) <= Shape::eps && dir[1] >= 0.0) return 0.0;
        if (fabs(p[2]-zMin()) <= Shape::eps && dir[2] <= 0.0) return 0.0;
        if (fabs(p[2]-zMax()) <= Shape::eps && dir[2] >= 0.0) return 0.0;
    }
    
    double d;
    double s = std::numeric_limits<double>::max();
    auto distX = [&](double l) { return p[0] + dir[0]*l; };
    auto distY = [&](double l) { return p[1] + dir[1]*l; };
    auto distZ = [&](double l) { return p[2] + dir[2]*l; };
    auto isValidDistX = [&](double l) { return distX(l) >= xMin() && distX(l) <= xMax(); };
    auto isValidDistY = [&](double l) { return distY(l) >= yMin() && distY(l) <= yMax(); };
    auto isValidDistZ = [&](double l) { return distZ(l) >= zMin() && distZ(l) <= zMax(); };

    if (dir[0] != 0.0){
        double dest = NAN;
        if (p[0]-xMin() < -Shape::eps && dir[0] > 0.0) dest = xMin(); // less than xMin
        else if (fabs(p[0]-xMin()) <= Shape::eps && dir[0] > 0.0) dest = xMax(); // on xMin
        else if (p[0]-xMin() > Shape::eps && p[0]-xMax() < -Shape::eps) dest = dir[0] > 0.0 ? xMax() : xMin(); // in between xMin and xMax
        else if (fabs(p[0]-xMax()) <= Shape::eps && dir[0] < 0.0) dest = xMin(); // on xMax
        else if (p[0]-xMax() > Shape::eps && dir[0] < 0.0) dest = xMax(); // greater than xMax

        if (std::isnan(dest)) return NAN;
        else{
            d = (dest-p[0])/dir[0];
            if (isValidDistY(d) && isValidDistZ(d)){ s = (d < s) ? d : s; }
        }
    }

    if (dir[1] != 0.0){
        double dest = NAN;
        if (p[1]-yMin() < -Shape::eps && dir[1] > 0.0) dest = yMin(); // less than yMin
        else if (fabs(p[1]-yMin()) <= Shape::eps && dir[1] > 0.0) dest = yMax(); // on yMin
        else if (p[1]-yMin() > Shape::eps && p[1]-xMax() < -Shape::eps) dest = dir[1] > 0.0 ? yMax() : yMin(); // in between yMin and yMax
        else if (fabs(p[1]-yMax()) <= Shape::eps && dir[1] < 0.0) dest = yMin(); // on yMax
        else if (p[1]-yMax() > Shape::eps && dir[1] < 0.0) dest = yMax(); // greater than yMax

        if (std::isnan(dest)) return NAN;
        else{
            d = (dest-p[1])/dir[1];
            if (isValidDistZ(d) && isValidDistX(d)){ s = (d < s) ? d : s; }
        }
    }

    if (dir[2] != 0.0){
        double dest = NAN;
        if (p[2]-zMin() < -Shape::eps && dir[2] > 0.0) dest = zMin(); // less than zMin
        else if (fabs(p[2]-zMin()) <= Shape::eps && dir[2] > 0.0) dest = zMax(); // on zMin
        else if (p[2]-zMin() > Shape::eps && p[2]-zMax() < -Shape::eps) dest = dir[2] > 0.0 ? zMax() : zMin(); // in between zMin and zMax
        else if (fabs(p[2]-zMax()) <= Shape::eps && dir[2] < 0.0) dest = zMin(); // on zMax
        else if (p[2]-zMax() > Shape::eps && dir[2] < 0.0) dest = zMax(); // greater than zMax

        if (std::isnan(dest)) return NAN;
        else{
            d = (dest-p[2])/dir[2];
            if (isValidDistX(d) && isValidDistY(d)){ s = (d < s) ? d : s; }
        }
    }
    return s == std::numeric_limits<double>::max() ? NAN : s;
}

Eigen::Vector3d Box::centroid() const noexcept{
    double x = (xMin() + xMax())/2;
    double y = (yMin() + yMax())/2;
    double z = (zMin() + zMax())/2;
    return Eigen::Vector3d(x, y, z);
}

UnitVector3d Box::normal(const Eigen::Vector3d& pos) const{
    if (!surfaceContains(pos)) throw std::invalid_argument("ERROR: Eigen::Vector3d is not on the surface.");
    if (fabs(pos[0] - xMin()) < Shape::eps) return UnitVector3d(-1, 0, 0);
    if (fabs(pos[0] - xMax()) < Shape::eps) return UnitVector3d(1, 0, 0);
    if (fabs(pos[1] - yMin()) < Shape::eps) return UnitVector3d(0, -1, 0);
    if (fabs(pos[1] - yMax()) < Shape::eps) return UnitVector3d(0, 1, 0);
    if (fabs(pos[2] - zMin()) < Shape::eps) return UnitVector3d(0, 0, -1);
    return UnitVector3d(0, 0, 1);
}

std::ostream& Box::print(std::ostream& os) const noexcept{
    os << "[" << xMin() << ", " << xMax() << "] * ";
    os << "[" << yMin() << ", " << yMax() << "] * ";
    os << "[" << zMin() << ", " << zMax() << "]";
    return os;
}