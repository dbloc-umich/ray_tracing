/**
 * This property class is used to describe parametric surfaces, which are 3D surfaces mapped from R2.
*/

#ifndef PARAMETRIC_SURFACE_PROPERTIES_H
#define PARAMETRIC_SURFACE_PROPERTIES_H

#include "FunctionShapeProperties.h"

class ParametricSurface;
class ParametricSurfaceProperties : public FunctionShapeProperties<std::function<double(const Eigen::Vector2d&)>>{
    public:
    friend class ParametricSurface;
    using Hfunction_type = std::function<Eigen::Matrix2d(const Eigen::Vector2d&)>;
    explicit ParametricSurfaceProperties(const function_type& x, const function_type& y, const function_type& z,
                                         double uMin, double uMax, double vMin, double vMax,
                                         unsigned short uSymmetry=1, unsigned short vSymmetry=1,
                                         const function_type& xu=nullptr, const function_type& yu=nullptr, const function_type& zu=nullptr,
                                         const function_type& xv=nullptr, const function_type& yv=nullptr, const function_type& zv=nullptr,
                                         const Hfunction_type& Hx=nullptr, const Hfunction_type& Hy=nullptr, const Hfunction_type& Hz=nullptr);

    Eigen::Vector3d r(const Eigen::Vector2d& u) const{ return {_x(u), _y(u), _z(u)}; }
    Eigen::Vector3d ru(const Eigen::Vector2d& u) const{ return {_xu(u), _yu(u), _zu(u)}; }
    Eigen::Vector3d rv(const Eigen::Vector2d& u) const{ return {_xv(u), _yv(u), _zv(u)}; }

    Eigen::Vector2d gradx(const Eigen::Vector2d& u) const{ return {_xu(u), _xv(u)}; }
    Eigen::Vector2d grady(const Eigen::Vector2d& u) const{ return {_yu(u), _yv(u)}; }
    Eigen::Vector2d gradz(const Eigen::Vector2d& u) const{ return {_zu(u), _zv(u)}; }

    protected:
    function_type _x, _y, _z;
    double _u0, _u1, _v0, _v1;
    unsigned short _uSym, _vSym; // number of intervals is u or v symmetric over
    function_type _xu, _yu, _zu;
    function_type _xv, _yv, _zv;
    Hfunction_type _Hx, _Hy, _Hz;
    void checkClosedSurface() override;
    void computeExtrema() override;
    void computeSurfaceArea() override;
    void computeVolume() override;

    bool isPeriodic(Eigen::Index i, Eigen::Index j, const Eigen::VectorXd&);
    bool isDegenerate(Eigen::Index i, Eigen::Index j, const Eigen::VectorXd&);
};

#endif