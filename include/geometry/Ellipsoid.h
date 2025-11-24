#ifndef ELLIPSOID_H
#define ELLIPSOID_H
#include "Shape.h"

class Ellipsoid: public Shape{  
    public:
    // Axis-aligned ellipsoids
    Ellipsoid(const Eigen::Vector3d& pt, double a, double b, double c, std::shared_ptr<Material> mat=nullptr);
    Ellipsoid(double x, double y, double z, double a, double b, double c, std::shared_ptr<Material> mat=nullptr):
        Ellipsoid(Eigen::Vector3d(x, y, z), a, b, c, mat) {}
    // General ellipsoids
    Ellipsoid(const Eigen::Vector3d& pt, const Eigen::Matrix3d& M, std::shared_ptr<Material> mat=nullptr);
    Ellipsoid(double x, double y, double z, const Eigen::Matrix3d& M, std::shared_ptr<Material> mat=nullptr):
        Ellipsoid(Eigen::Vector3d(x, y, z), M, mat) {}

    Eigen::Vector3d center() const noexcept{ return _center; }
    double semiAxis(Eigen::Index i) const noexcept;
    UnitVector3d principalAxis(Eigen::Index i) const noexcept;

    void setCenter(const Eigen::Vector3d& point){ _center = point; }
    void setSemiAxis(Eigen::Index i, double L);

    double xMin() const noexcept override;
    double xMax() const noexcept override;
    double yMin() const noexcept override;
    double yMax() const noexcept override;
    double zMin() const noexcept override;
    double zMax() const noexcept override;

    double surfaceArea() const noexcept override;
    double volume() const noexcept override;

    bool surfaceContains(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Eigen::Vector3d& p) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;
    double distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept override;

    UnitVector3d normal(const Eigen::Vector3d& pos) const override;

    protected:
    Eigen::Vector3d _center;
    Eigen::MatrixXd _M;
    Eigen::MatrixXd _Minv; // inverse of M
    Eigen::MatrixXd _axes; // eigenvectors of _M, principal axis directions
    Eigen::Array3d _length; // inverse squares of eigenvalues of _M
    std::ostream& print(std::ostream& os) const noexcept override;

    double xdotMy(const Eigen::Vector3d&, const Eigen::Vector3d&) const noexcept;
    bool intersects(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1) const noexcept;
    bool intersects(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) const noexcept;
};

#endif // ELLIPSOID_H