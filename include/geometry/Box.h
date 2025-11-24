#ifndef BOX_H
#define BOX_H
#include "Shape.h"

class Box: public Shape{
    public:
    Box(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper, std::shared_ptr<Material> mat = nullptr);
    explicit Box(double x0=0.0, double y0=0.0, double z0=0.0, double x1=1.0, double y1=1.0, double z1=1.0,
                 std::shared_ptr<Material> mat = nullptr);

    Eigen::Vector3d lowerVertex() const noexcept{ return _lower; }
    Eigen::Vector3d upperVertex() const noexcept{ return _upper; }

    void setLowerVertex(const Eigen::Vector3d& point);
    void setUpperVertex(const Eigen::Vector3d& point);
    void setVertices(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper);

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

    bool surfaceContains(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Shape& other) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;
    double distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept override;

    UnitVector3d normal(const Eigen::Vector3d& pos) const override;

    protected:
    Eigen::Vector3d _lower, _upper;
    std::ostream& print(std::ostream& os) const noexcept override;
};
#endif // BOX_H