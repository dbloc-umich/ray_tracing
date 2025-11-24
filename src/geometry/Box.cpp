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
    bool onX = std::abs(p[0]-xMin()) <= Shape::eps || std::abs(p[0]-xMax()) <= Shape::eps;
    bool onY = std::abs(p[1]-yMin()) <= Shape::eps || std::abs(p[1]-yMax()) <= Shape::eps;
    bool onZ = std::abs(p[2]-zMin()) <= Shape::eps || std::abs(p[2]-zMax()) <= Shape::eps;
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
    }
    return other.overlaps(*this); // See the implementation files for other Shapes
}

double Box::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{
    bool isOnSurface = surfaceContains(p);
    if (isOnSurface){ // point moving on a surface
        if (std::abs(p[0]-xMin()) <= Shape::eps && dir[0] <= 0.0) return 0.0;
        if (std::abs(p[0]-xMax()) <= Shape::eps && dir[0] >= 0.0) return 0.0;
        if (std::abs(p[1]-yMin()) <= Shape::eps && dir[1] <= 0.0) return 0.0;
        if (std::abs(p[1]-yMax()) <= Shape::eps && dir[1] >= 0.0) return 0.0;
        if (std::abs(p[2]-zMin()) <= Shape::eps && dir[2] <= 0.0) return 0.0;
        if (std::abs(p[2]-zMax()) <= Shape::eps && dir[2] >= 0.0) return 0.0;
    }
    
    Eigen::Array<double, 6, 1> s;
    s[0] = (xMin()-p[0])/dir[0] > 0 ? (xMin()-p[0])/dir[0] : std::numeric_limits<double>::quiet_NaN(); // distance to x = xMin plane
    s[1] = (xMax()-p[0])/dir[0] > 0 ? (xMax()-p[0])/dir[0] : std::numeric_limits<double>::quiet_NaN(); // distance to x = xMax plane
    s[2] = (yMin()-p[1])/dir[1] > 0 ? (yMin()-p[1])/dir[1] : std::numeric_limits<double>::quiet_NaN(); // distance to y = yMin plane
    s[3] = (yMax()-p[1])/dir[1] > 0 ? (yMax()-p[1])/dir[1] : std::numeric_limits<double>::quiet_NaN(); // distance to y = yMax plane
    s[4] = (zMin()-p[2])/dir[2] > 0 ? (zMin()-p[2])/dir[2] : std::numeric_limits<double>::quiet_NaN(); // distance to z = zMin plane
    s[5] = (zMax()-p[2])/dir[2] > 0 ? (yMax()-p[2])/dir[2] : std::numeric_limits<double>::quiet_NaN(); // distance to z = zMax plane

    // Checking the positive distances to see whether they are still valid
    for (Eigen::Index i = 0; i < 6; i++){
        if (!std::isnan(s[i])){
            Eigen::Vector3d dest = p + dir.value()*s[i];
            if (i <= 1){ // Checking if y and z are in bound
                if (dest[1] < yMin() || dest[1] > yMax() || dest[2] < zMin() || dest[2] > zMax()){
                    s[i] = std::numeric_limits<double>::quiet_NaN();
                }
            } else if (i <= 3){ // Checking if x and z are in bound
                if (dest[0] < xMin() || dest[0] > xMax() || dest[2] < zMin() || dest[2] > zMax()){
                    s[i] = std::numeric_limits<double>::quiet_NaN();
                }
            } else{ // Checking if x and y are in bound
                if (dest[0] < xMin() || dest[0] > xMax() || dest[1] < yMin() || dest[1] > yMax()){
                    s[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }

    if (s.isNaN().all()) return s[0]; // return a NaN
    double dist = std::numeric_limits<double>::max();
    for (double d: s){
        if (dist > d) dist = d;
    }
    return dist;
}

UnitVector3d Box::normal(const Eigen::Vector3d& pos) const{
    // Count how many faces the point lies on
    short s = std::abs(pos[0] - xMin()) < Shape::eps;
    s += std::abs(pos[0] - xMax()) < Shape::eps;
    s += std::abs(pos[1] - yMin()) < Shape::eps;
    s += std::abs(pos[1] - yMax()) < Shape::eps;
    s += std::abs(pos[2] - zMin()) < Shape::eps;
    s += std::abs(pos[2] - zMax()) < Shape::eps;

    if (s == 0) throw std::invalid_argument("ERROR: Eigen::Vector3d is not on the surface.");
    if (s > 1) return UnitVector3d(nullptr); // lies on multiple faces

    if (std::abs(pos[0] - xMin()) < Shape::eps) return UnitVector3d(-1, 0, 0);
    if (std::abs(pos[0] - xMax()) < Shape::eps) return UnitVector3d(1, 0, 0);
    if (std::abs(pos[1] - yMin()) < Shape::eps) return UnitVector3d(0, -1, 0);
    if (std::abs(pos[1] - yMax()) < Shape::eps) return UnitVector3d(0, 1, 0);
    if (std::abs(pos[2] - zMin()) < Shape::eps) return UnitVector3d(0, 0, -1);
    return UnitVector3d(0, 0, 1);
}

std::ostream& Box::print(std::ostream& os) const noexcept{
    os << "[" << xMin() << ", " << xMax() << "] * ";
    os << "[" << yMin() << ", " << yMax() << "] * ";
    os << "[" << zMin() << ", " << zMax() << "]";
    return os;
}