#ifndef RAY_H
#define RAY_H

#include "UnitVector3d.h"

class Shape;
class Ray{
    public:
    Ray(Eigen::Vector3d p, UnitVector3d dir, double I=1.0, double lambda=325e-9, Shape* host=nullptr);

    Eigen::Vector3d position() const noexcept{ return _p; }
    void setPoistion(const Eigen::Vector3d& p) noexcept{ _p = p; }

    UnitVector3d direction() const noexcept{ return _dir; }
    void setDirection(const UnitVector3d& dir) noexcept{ _dir = dir; }

    double intensity() const noexcept{ return _I; }
    void setIntensity(double I);

    double wavelength() const noexcept{ return _lambda; }
    double frequency() const noexcept;
    void setWavelength(double lambda);

    Shape* host() const noexcept{ return _host; }
    void setHost(Shape* host);

    protected:
    Eigen::Vector3d _p;
    UnitVector3d _dir;
    double _I;
    double _lambda; // vacuum wavelength, [m]
    Shape* _host;
};

#endif // RAY_H