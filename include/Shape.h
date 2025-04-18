#ifndef SHAPE_H
#define SHAPE_H
#include <iostream>
#include <limits>

class Direction;
class Point;

class Shape{
    public:
    explicit Shape(double Sigma_t=0.0, double refrac=1.0);
    virtual ~Shape() = default;

    virtual double Sigma_t() const noexcept{ return _Sigma_t; }
    void setSigma_t(double Sigma_t);

    virtual double refractive() const noexcept{ return _refrac; }
    void setRefractive(double refrac);

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
    virtual bool contentsOverlap(const Shape& other) const noexcept{ return this->overlaps(other); };
    /***
     * if the Point is inside the Shape, the smallest (positive) distance to the surface(s)
     * if the Point is on the surface, 0.0 if it leaves and the smallest distance to the surface if it enters
     * if the Point is outside the Shape, NAN if it never enters the Shape, else the smallest distance to the surface
    ***/
    virtual double distanceToSurface(const Point& p, const Direction& dir) const noexcept = 0;

    virtual Direction normal(const Point& pos) const = 0; // outward unit normal vector
    friend std::ostream& operator<<(std::ostream& os, const Shape& shape){ return shape.print(os); }

    protected:
    virtual std::ostream& print(std::ostream& os) const noexcept = 0;
    static constexpr double eps = 1e-9;
    double _Sigma_t, _refrac;
};
#endif // SHAPE_H