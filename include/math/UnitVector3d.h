#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include "Eigen/Dense"

class UnitVector3d{
    public:
    explicit UnitVector3d(std::nullptr_t);
    explicit UnitVector3d(const Eigen::Vector3d& v);
    UnitVector3d(double x, double y, double z): UnitVector3d(Eigen::Vector3d(x, y, z)) {}

    Eigen::Vector3d& value() noexcept{ return _v; }
    const Eigen::Vector3d& value() const noexcept{ return _v; }
    double operator[](std::size_t i) const noexcept{ return _v[i]; }
    explicit operator Eigen::Vector3d() const noexcept{ return _v; }
    explicit operator bool() const noexcept;

    bool isOrthogonal(const UnitVector3d& other) const noexcept;
    bool isParallel(const UnitVector3d& other) const noexcept;

    UnitVector3d reflect(const UnitVector3d& normal);
    UnitVector3d refract(const UnitVector3d& normal, double n1, double n2);

    friend std::ostream& operator<<(std::ostream& os, const UnitVector3d& v);

    protected:
    Eigen::Vector3d _v;
};

#endif // VECTOR_H