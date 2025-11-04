#ifndef SHAPE_H
#define SHAPE_H

#include <memory>
#include <vector>
#include "UnitVector3d.h"

class Material;
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

    virtual bool surfaceContains(const Eigen::Vector3d& p) const noexcept = 0;
    virtual bool encloses(const Eigen::Vector3d& p) const noexcept = 0; // If this->surfaceContains(p) is true, then this->encloses(p) is false;
    virtual bool encloses(const Shape& other) const noexcept = 0;
    virtual bool overlaps(const Shape& other) const noexcept = 0;
    /***
     * if the Point is inside the Shape, the smallest (positive) distance to the surface(s)
     * if the Point is on the surface, 0.0 if it leaves and the smallest distance to the surface if it enters
     * if the Point is outside the Shape, NAN if it never enters the Shape, else the smallest distance to the surface
    ***/
    virtual double distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept = 0;

    virtual Eigen::Vector3d centroid() const noexcept = 0;
    virtual UnitVector3d normal(const Eigen::Vector3d& pos) const = 0; // outward unit normal vector
    friend std::ostream& operator<<(std::ostream& os, const Shape& shape);

    bool hasProperty(const Prop&) const noexcept;
    double computeProperty(const Prop&, const std::vector<double>& vars={}) const;

    protected:
    std::shared_ptr<Material> _mat;
    static constexpr double eps = 1e-9;
    virtual std::ostream& print(std::ostream& os) const noexcept = 0;
};
#endif // SHAPE_H