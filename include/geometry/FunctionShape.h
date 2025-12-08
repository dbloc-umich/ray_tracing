#ifndef FUNCTION_SHAPE_H
#define FUNCTION_SHAPE_H

#include "FunctionShapeProperties.h"
#include "Shape.h"
#include <Eigen/SVD>

template<typename T>
class FunctionShape : public Shape{
    public:
    using property_ptr = std::shared_ptr<T>;
    explicit FunctionShape(property_ptr prop,
                           const Eigen::Matrix3d& M = Eigen::Matrix3d::Identity(),
                           const Eigen::Vector3d& dr = Eigen::Vector3d::Zero(),
                           std::shared_ptr<Material> mat=nullptr):
        Shape(mat),
        _prop(prop),
        _M(M),
        _only_scaling(M.isApprox(M.diagonal().asDiagonal().toDenseMatrix())),
        _dr(dr)
    {
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(M);
        _sigma = svd.singularValues();
        if (_sigma[0] < Shape::eps || _sigma[1] < Shape::eps || _sigma[2] < Shape::eps){
            throw std::invalid_argument("ERROR: Singular matrix detected.");
        }
    };
    
    virtual ~FunctionShape() = default;
    double volume() const noexcept override{
        if (_prop->isClosedSurface()) return _prop->volume() * std::abs(_sigma[0]*_sigma[1]*_sigma[2]);
        return _prop->volume(); // should be NAN if the surface is open
    }

    protected:
    property_ptr _prop; // property pointer to a base Shape
    Eigen::Matrix3d _M; // linear transformation matrix
    bool _only_scaling; // whether _M is a scaling-only matrix
    Eigen::Vector3d _dr; // translation vector
    Eigen::Vector3d _sigma; // singular values of _M

    Eigen::ArrayXd _lower, _upper; // lower and upper vertex of the bounding box
    virtual void computeExtrema() = 0; // only needed if the extreme points cannot be easily computed from _prop
    double scaledXMin() const noexcept{ return _prop->xMin() * _M(0,0) + _dr[0]; }
    double scaledXMax() const noexcept{ return _prop->xMax() * _M(0,0) + _dr[0]; }
    double scaledYMin() const noexcept{ return _prop->yMin() * _M(1,1) + _dr[1]; }
    double scaledYMax() const noexcept{ return _prop->yMax() * _M(1,1) + _dr[1]; }
    double scaledZMin() const noexcept{ return _prop->zMin() * _M(2,2) + _dr[2]; }
    double scaledZMax() const noexcept{ return _prop->zMax() * _M(2,2) + _dr[2]; }
    
    double _surfaceArea;
    virtual void computeSurfaceArea() = 0; // only needed if the surface areas cannot be easily computed from _prop
    double uniformlyScaledSurfaceArea() const noexcept{ return _prop->surfaceArea() * (_sigma[0]*_sigma[0]); }
};

#endif