#ifndef RAY_H
#define RAY_H

#include "UnitVector.h"

class Shape;
class Ray{
    public:
    Ray(Eigen::Vector3d p, UnitVector3d dir, double I, double lambda, double A, Shape* host=nullptr);

    Eigen::Vector3d position() const noexcept{ return _p; }
    void setPoistion(const Eigen::Vector3d& p) noexcept{ _p = p; }

    UnitVector3d direction() const noexcept{ return _dir; }
    void setDirection(const UnitVector3d& dir) noexcept{ _dir = dir; }

    double intensity() const noexcept{ return _I; }
    void setIntensity(double I);

    double wavelength() const noexcept{ return _lambda; }
    double frequency() const noexcept; // angular frequency in radians/s
    void setWavelength(double lambda);

    double area() const noexcept{ return _A; }
    double power() const noexcept{ return _I*_A; }
    void setArea(double A);

    Shape* host() const noexcept{ return _host; }
    void setHost(Shape* host);

    protected:
    Eigen::Vector3d _p;
    UnitVector3d _dir;
    double _I; // intensity [W/m2]
    double _lambda; // vacuum wavelength, [m]
    double _A; // cross-section area, [m2]
    Shape* _host;
};

#endif // RAY_H