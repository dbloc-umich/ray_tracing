#ifndef SHAPE_H
#define SHAPE_H

#include <memory>
#include <vector>

class Direction;
class Material;
class Point;
enum class Prop;

class Shape{
    public:
    explicit Shape(std::shared_ptr<Material> mat=nullptr);
    virtual ~Shape() = default;

    virtual double xMin() const noexcept = 0;
    virtual double xMax() const noexcept = 0;
    virtual double yMin() const noexcept = 0;
    virtual double yMax() const noexcept = 0;
    virtual double zMin() const noexcept = 0;
    virtual double zMax() const noexcept = 0;

    virtual double surfaceArea() const noexcept = 0;
    virtual double volume() const noexcept = 0;

    virtual bool surfaceContains(const Point& p) const noexcept = 0;
    virtual bool encloses(const Point& p) const noexcept = 0; // If this->surfaceContains(p) is true, then this->encloses(p) is false;
    virtual bool encloses(const Shape& other) const noexcept = 0;
    virtual bool overlaps(const Shape& other) const noexcept = 0;
    /***
     * if the Point is inside the Shape, the smallest (positive) distance to the surface(s)
     * if the Point is on the surface, 0.0 if it leaves and the smallest distance to the surface if it enters
     * if the Point is outside the Shape, NAN if it never enters the Shape, else the smallest distance to the surface
    ***/
    virtual double distanceToSurface(const Point& p, const Direction& dir) const noexcept = 0;

    virtual Point centroid() const noexcept = 0;
    virtual Direction normal(const Point& pos) const = 0; // outward unit normal vector
    double getProp(const Prop&, const std::vector<double>& vars={}) const;
    friend std::ostream& operator<<(std::ostream& os, const Shape& shape);

    protected:
    std::shared_ptr<Material> _mat;
    static constexpr double eps = 1e-9;
    virtual std::ostream& print(std::ostream& os) const noexcept = 0;
};
#endif // SHAPE_H