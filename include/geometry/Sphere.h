#ifndef SPHERE_H
#define SPHERE_H
#include "Shape.h"

class Sphere: public Shape{  
    public:
    explicit Sphere(const Eigen::Vector3d& pt, double R=1.0, std::shared_ptr<Material> mat=nullptr);
    explicit Sphere(double x=0.0, double y=0.0, double z=0.0, double R=1.0,
                    std::shared_ptr<Material> mat = nullptr);

    Eigen::Vector3d origin() const noexcept{ return _origin; }
    double radius() const noexcept{ return _radius; }

    void setOrigin(const Eigen::Vector3d& point){ _origin = point; }
    void setRadius(double R);

    double xMin() const noexcept override{ return _origin.x() - _radius; }
    double xMax() const noexcept override{ return _origin.x() + _radius; }
    double yMin() const noexcept override{ return _origin.y() - _radius; }
    double yMax() const noexcept override{ return _origin.y() + _radius; }
    double zMin() const noexcept override{ return _origin.z() - _radius; }
    double zMax() const noexcept override{ return _origin.z() + _radius; }

    double surfaceArea() const noexcept override;
    double volume() const noexcept override;

    bool surfaceContains(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Shape& other) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;
    double distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept override;

    Eigen::Vector3d centroid() const noexcept override{ return _origin; }
    UnitVector3d normal(const Eigen::Vector3d& pos) const override;

    protected:
    Eigen::Vector3d _origin;
    double _radius;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif // SPHERE_H