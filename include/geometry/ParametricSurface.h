#ifndef PARAMETRIC_FUNCTION_SHAPE_H
#define PARAMETRIC_FUNCTION_SHAPE_H

#include "FunctionShape.h"
#include "ParametricSurfaceProperties.h"
#include <functional>

class ParametricSurface : public FunctionShape<ParametricSurfaceProperties>{
    public:
    explicit ParametricSurface(property_ptr prop,
                               const Eigen::Matrix3d& M = Eigen::Matrix3d::Identity(),
                               const Eigen::Vector3d& dr = Eigen::Vector3d::Zero(),
                               std::shared_ptr<Material> mat=nullptr);

    double xMin() const noexcept override{ return _only_scaling ? scaledXMin() : _lower[0]; };
    double xMax() const noexcept override{ return _only_scaling ? scaledXMax() : _upper[0]; };
    double yMin() const noexcept override{ return _only_scaling ? scaledYMin() : _lower[1]; };
    double yMax() const noexcept override{ return _only_scaling ? scaledYMax() : _upper[1]; };
    double zMin() const noexcept override{ return _only_scaling ? scaledZMin() : _lower[2]; };
    double zMax() const noexcept override{ return _only_scaling ? scaledZMax() : _upper[2]; };
    double surfaceArea() const noexcept override{ return _sigma[0]/_sigma[2] - 1 <= Shape::eps ? uniformlyScaledSurfaceArea() : _surfaceArea; }

    bool surfaceContains(const Eigen::Vector3d& p) const noexcept override;
    bool encloses(const Eigen::Vector3d& p) const noexcept override;
    bool overlaps(const Shape& other) const noexcept override;

    double distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept override;
    UnitVector3d normal(const Eigen::Vector3d& pos) const override; // outward unit normal vector

    Eigen::Vector3d r(const Eigen::Vector2d& u) const{ return _M*_prop->r(u)+_dr; }
    Eigen::Vector3d ru(const Eigen::Vector2d& u) const{ return _M*_prop->ru(u); }
    Eigen::Vector3d rv(const Eigen::Vector2d& u) const{ return _M*_prop->rv(u); }

    protected:
    void computeExtrema() override;
    void computeSurfaceArea() override;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif