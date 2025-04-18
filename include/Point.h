#ifndef POINT_H
#define POINT_H
#include <iostream>

class Vector;
class Point{
    public:
    explicit Point(double x=0.0, double y=0.0, double z=0.0): _x(x), _y(y), _z(z) {}

    double x() const noexcept{ return _x; }
    double y() const noexcept{ return _y; }
    double z() const noexcept{ return _z; }

    void setX(double x) noexcept{ _x = x; }
    void setY(double y) noexcept{ _y = y; }
    void setZ(double z) noexcept{ _z = z; }

    double dist(const Point& other) const noexcept;
    bool operator==(const Point& other) const noexcept;
    bool operator!=(const Point& other) const noexcept{ return !(*this==other); }

    Point& advance(const Vector& v, const double s=1);
    
    friend std::ostream& operator<<(std::ostream& os, const Point& p);

    private:
    double _x, _y, _z;
};
#endif // POINT_H