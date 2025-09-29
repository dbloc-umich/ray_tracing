#ifndef SPHERE_H
#define SPHERE_H
#include "Shape.h"
#include "Point.h"

class Box;
class Sphere: public Shape{  
    public:
    explicit Sphere(const Point& pt, double R=1.0, double Sigma_t=0.0, double refrac=1.0);
    explicit Sphere(double x=0.0, double y=0.0, double z=0.0, double R=1.0, double Sigma_t=0.0, double refrac=1.0);

    Point origin() const noexcept{ return _origin; }
    double radius() const noexcept{ return _radius; }

    void setOrigin(const Point& point){ _origin = point; }
    void setRadius(double R);

    double xMin() const noexcept override{ return _origin.x() - _radius; }
    double xMax() const noexcept override{ return _origin.x() + _radius; }
    double yMin() const noexcept override{ return _origin.y() - _radius; }
    double yMax() const noexcept override{ return _origin.y() + _radius; }
    double zMin() const noexcept override{ return _origin.z() - _radius; }
    double zMax() const noexcept override{ return _origin.z() + _radius; }

    double surfaceArea() const noexcept override;
    double volume() const noexcept override;

    bool surfaceContains(const Point& p) const noexcept override;
    bool encloses(const Point& p) const noexcept override;
    bool encloses(const Shape& other) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;
    double distanceToSurface(const Point& p, const Direction& dir) const noexcept override;

    Point centroid() const noexcept override{ return _origin; }
    Direction normal(const Point& pos) const override;

    protected:
    Point _origin;
    double _radius;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif // SPHERE_H