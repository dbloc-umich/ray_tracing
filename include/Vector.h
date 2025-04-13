#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>

class Point;
class Vector{
    public:
    Vector(double dx, double dy, double dz);
    Vector(const Point& p1, const Point& p2);

    double dx() const noexcept{ return _dx; }
    double dy() const noexcept{ return _dy; }
    double dz() const noexcept{ return _dz; }

    void setX(double dx) noexcept{ _dx = dx; }
    void setY(double dy) noexcept{ _dy = dy; }
    void setZ(double dz) noexcept{ _dz = dz; }

    Vector operator+() const noexcept;
    Vector operator+(const Vector& other) const noexcept;
    Vector& operator+=(const Vector& other) noexcept;

    Vector operator-() const noexcept;
    Vector operator-(const Vector& other) const noexcept;
    Vector& operator-=(const Vector& other) noexcept;

    Vector operator*(double d) const noexcept;
    Vector& operator*=(double d) noexcept;

    Vector operator/(double d) const;
    Vector& operator/=(double d);

    bool operator==(const Vector& other) const noexcept;
    bool operator!=(const Vector& other) const noexcept{ return !(*this==other); }
    explicit operator bool() const noexcept;

    double dot(const Vector& other) const noexcept;
    virtual double norm() const noexcept;
    Vector cross(const Vector& other) const noexcept;

    bool isOrthogonal(const Vector& other) const noexcept;
    bool isParallel(const Vector& other) const noexcept;

    friend std::ostream& operator<<(std::ostream& os, const Vector& p);

    protected:
    double _dx, _dy, _dz;
};

inline Vector operator*(double d, const Vector& v) noexcept{ return v*d; }
inline Vector operator/(double d, const Vector& v){ return v/d; }

#endif // VECTOR_H