#ifndef BOX_H
#define BOX_H
#include "Shape.h"
#include "Point.h"

class Sphere;
class Box: public Shape{
    public:
    Box(const Point& lower, const Point& upper, double Sigma_t=0.0, double refrac=1.0);
    explicit Box(double x0=0.0, double y0=0.0, double z0=0.0, double x1=1.0, double y1=1.0, double z1=1.0, double Sigma_t=0.0, double refrac=1.0);

    Point lowerVertex() const noexcept{ return _lower; }
    Point upperVertex() const noexcept{ return _upper; }

    void setLowerVertex(const Point& point);
    void setUpperVertex(const Point& point);
    void setVertices(const Point& lower, const Point& upper);

    double xMin() const noexcept override{ return _lower.x(); }
    double xMax() const noexcept override{ return _upper.x(); }
    double yMin() const noexcept override{ return _lower.y(); }
    double yMax() const noexcept override{ return _upper.y(); }
    double zMin() const noexcept override{ return _lower.z(); }
    double zMax() const noexcept override{ return _upper.z(); }

    double length() const noexcept{ return _upper.x()-_lower.x(); }
    double width() const noexcept{ return _upper.y()-_lower.y(); }
    double height() const noexcept{ return _upper.z()-_lower.z(); }

    double surfaceArea() const noexcept override{ return 2*(length()*width() + length()*height() + width()*height()); }
    double volume() const noexcept override{ return length()*width()*height(); }

    bool surfaceContains(const Point& p) const noexcept override;
    bool encloses(const Point& p) const noexcept override;
    bool encloses(const Shape& other) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;
    double distanceToSurface(const Point& p, const Direction& dir) const noexcept override;

    Point centroid() const noexcept override;
    Direction normal(const Point& pos) const override;

    protected:
    Point _lower, _upper;
    std::ostream& print(std::ostream& os) const noexcept override;
};
#endif // BOX_H