#ifndef FV_SPHERICAL_MESH_H
#define FV_SPHERICAL_MESH_H

#include "FVSpatialMesh.h"

class Sphere;
class FVSphericalMesh: public FVSpatialMesh{
    public:
    FVSphericalMesh(const Eigen::ArrayXd& r, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& phi);
    FVSphericalMesh(const std::array<Eigen::ArrayXd, 3>& axes);

    bool isOrthogonal() const noexcept override{ return true; };
    Eigen::Vector3d cartesian(double, double, double) const noexcept override;
    Eigen::Vector3d coordinates(double, double, double) const noexcept override;
    Eigen::Vector3d basisVector(Eigen::Index varID, double r, double mu, double phi) const noexcept override;
    UnitVector3d normalVector(Eigen::Index varID, double r, double mu, double phi) const noexcept override{ return basisVector(varID, r, mu, phi); }
    double gradientFactor(Eigen::Index varID, double r, double mu, double phi) const noexcept override;
    double gradientDotN(double dudq, Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const override;
    double area(Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const noexcept override;
    double volume(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept override;
    double ghostLength(Eigen::Index surfID) const noexcept override;

    protected:
    Eigen::ArrayXd& _r, _mu, _phi; // alias for existing axes
};

#endif