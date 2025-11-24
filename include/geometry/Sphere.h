#ifndef SPHERE_H
#define SPHERE_H
#include "Shape.h"

class Sphere: public Shape{  
    public:
    explicit Sphere(const Eigen::Vector3d& pt, double R=1.0, std::shared_ptr<Material> mat=nullptr);
    explicit Sphere(double x=0.0, double y=0.0, double z=0.0, double R=1.0,
                    std::shared_ptr<Material> mat = nullptr);

    Eigen::Vector3d center() const noexcept{ return _center; }
    double radius() const noexcept{ return _radius; }

    void setCenter(const Eigen::Vector3d& point){ _center = point; }
    void setRadius(double R);

    double xMin() const noexcept override{ return _center.x() - _radius; }
    double xMax() const noexcept override{ return _center.x() + _radius; }
    double yMin() const noexcept override{ return _center.y() - _radius; }
    double yMax() const noexcept override{ return _center.y() + _radius; }
    double zMin() const noexcept override{ return _center.z() - _radius; }
    double zMax() const noexcept override{ return _center.z() + _radius; }

    double surfaceArea() const noexcept override;
    double volume() const noexcept override;

    bool surfaceContains(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Eigen::Vector3d& p) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;
    double distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept override;

    UnitVector3d normal(const Eigen::Vector3d& pos) const override;

    protected:
    Eigen::Vector3d _center;
    double _radius;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif // SPHERE_H