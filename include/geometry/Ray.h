#ifndef RAY_H
#define RAY_H

#include "Direction.h"
#include "Point.h"

class Shape;
class Ray{
    public:
    Ray(Point p, Direction dir, double I=1.0, double lambda=325e-9, Shape* host=nullptr);

    Point position() const noexcept{ return _p; }
    void setPoistion(const Point& p) noexcept{ _p = p; }

    Direction direction() const noexcept{ return _dir; }
    void setDirection(const Direction& dir) noexcept{ _dir = dir; }

    double intensity() const noexcept{ return _I; }
    void setIntensity(double I);

    double wavelength() const noexcept{ return _lambda; }
    void setWavelength(double lambda);

    Shape* host() const noexcept{ return _host; }
    void setHost(Shape* host);

    protected:
    Point _p;
    Direction _dir;
    double _I;
    double _lambda;
    Shape* _host;
};

#endif // RAY_H