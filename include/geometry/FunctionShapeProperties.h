/**
 * The motivation for this class is so that a FuntionShape object can have a shared pointer to this. Many of these properties may be expensive to compute, so having a
 * pointer to this property object will help reduce redundant calculations. For identical objects that differ only in their position, properties such as surface area and
 * volume do not change, while other only differ by the same displacement.
*/

#ifndef FUNCTION_SHAPE_PROPERTIES_H
#define FUNCTION_SHAPE_PROPERTIES_H
#include "Eigen/Dense"
#include <type_traits>

template<typename Callable,
         typename = std::enable_if_t<!std::is_reference<Callable>::value>>
class FunctionShapeProperties{
    public:
    using function_type = Callable;
    virtual ~FunctionShapeProperties() = default;

    double xMin() const noexcept{ return _xMin; }
    double xMax() const noexcept{ return _xMax; }
    double yMin() const noexcept{ return _yMin; }
    double yMax() const noexcept{ return _yMax; }
    double zMin() const noexcept{ return _zMin; }
    double zMax() const noexcept{ return _zMax; }
    double surfaceArea() const noexcept{ return _surfaceArea; }
    double volume() const noexcept{ return _volume; }
    bool isClosedSurface() const noexcept{ return _isClosedSurface; }

    protected:
    // These variables are computed upon construction
    bool _isClosedSurface;
    double _xMin;
    double _xMax;
    double _yMin;
    double _yMax;
    double _zMin;
    double _zMax;
    double _surfaceArea;
    double _volume;

    // Required helper functions
    virtual void checkClosedSurface() = 0;
    virtual void computeExtrema() = 0;
    virtual void computeSurfaceArea() = 0;
    virtual void computeVolume() = 0;
    static constexpr double eps = 1.0e-6;
};

#endif