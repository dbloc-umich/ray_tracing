#ifndef FV_ELLIPSOIDAL_MESH_H
#define FV_ELLIPSOIDAL_MESH_H

#include "FVSpatialMesh.h"

class Ellipsoid;
class FVEllipsoidalMesh: public FVSpatialMesh{
    public:
    FVEllipsoidalMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, double a, double b, double c);
    FVEllipsoidalMesh(const std::array<Eigen::ArrayXd, 3>& axes, double a, double b, double c);
    FVEllipsoidalMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, const Ellipsoid& ell);
    FVEllipsoidalMesh(const std::array<Eigen::ArrayXd, 3>& axes, const Ellipsoid& ell);

    bool isOrthogonal() const noexcept override{ return (_a == _b && _a == _c); };
    Eigen::Vector3d cartesian(double, double, double) const noexcept override;
    Eigen::Vector3d coordinates(double, double, double) const noexcept override;
    Eigen::Vector3d basisVector(Eigen::Index varID, double r, double mu, double phi) const noexcept override;
    UnitVector3d normalVector(Eigen::Index varID, double r, double mu, double phi) const noexcept override;
    double area(Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const noexcept override;
    double volume(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept override; 

    protected:
    double _a, _b, _c;
    Eigen::ArrayXd& _r, _mu, _phi; // alias for existing axes
};

#endif