#include "FVSphericalMesh.h"
#include "Constants.h"
#include "Sphere.h"

FVSphericalMesh::FVSphericalMesh(const Eigen::ArrayXd& r, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& phi):
    FVSpatialMesh(r, mu, phi),
    _r(_axes[0]),
    _mu(_axes[1]),
    _phi(_axes[2])
{
    if (_r[0] < 0) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_mu[0] < -1 || _mu[_mu.size()-1] > 1) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_phi[0] < -mconst::pi || _phi[_phi.size()-1] > mconst::pi) throw std::invalid_argument("ERROR: Invalid boundary values.");
}

FVSphericalMesh::FVSphericalMesh(const std::array<Eigen::ArrayXd, 3>& axes):
    FVSpatialMesh(axes),
    _r(_axes[0]),
    _mu(_axes[1]),
    _phi(_axes[2])    
{
    if (_r[0] < 0) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_mu[0] < -1 || _mu[_mu.size()-1] > 1) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_phi[0] < -mconst::pi || _phi[_phi.size()-1] > mconst::pi) throw std::invalid_argument("ERROR: Invalid boundary values.");
}

Eigen::Vector3d FVSphericalMesh::cartesian(double r, double mu, double phi) const noexcept{
    return {r*mu, r*std::sqrt(1-mu*mu)*std::cos(phi), r*std::sqrt(1-mu*mu)*std::sin(phi)};
}

Eigen::Vector3d FVSphericalMesh::coordinates(double x, double y, double z) const noexcept{
    double r = std::sqrt(x*x + y*y + z*z);
    return {r, x/r, std::atan2(z,y)};
}

Eigen::Vector3d FVSphericalMesh::basisVector(Eigen::Index varID, double r, double mu, double phi) const noexcept{
    // varID = 0, 1, 2 referring to which basis vector to return
    if (varID == 0) return {mu, std::sqrt(1-mu*mu)*std::cos(phi), std::sqrt(1-mu*mu)*std::sin(phi)};
    if (varID == 1) return {r, -r*mu/std::sqrt(1-mu*mu)*std::cos(phi), -r*mu/std::sqrt(1-mu*mu)*std::sin(phi)};
    return {0, -r*std::sqrt(1-mu*mu)*std::sin(phi), r*std::sqrt(1-mu*mu)*std::cos(phi)};
}

double FVSphericalMesh::area(Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const noexcept{
    if (surfID <= 1){
        // surfaces where r is constant
        double r = surfID == 0 ? _r[i] : _r[i+1];
        if (r == 0.0) return 0;
        return r*r * (_mu[j+1]-_mu[j]) * (_phi[k+1]-_phi[k]);
    } else if (surfID <= 3){
        // surfaces where mu is constant
        double mu = surfID == 2 ? _mu[j] : _mu[j+1];
        if (mu == 1.0 || mu == -1.0) return 0;
        return (_r[i+1]*_r[i+1]-_r[i]*_r[i])/2 * std::sqrt(1-mu*mu) * (_phi[k+1]-_phi[k]);
    }
    // surfaces where phi is constant
    return (_r[i+1]*_r[i+1]-_r[i]*_r[i])/2 * (std::asin(_mu[j+1])-std::asin(_mu[j]));
}

double FVSphericalMesh::volume(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    double rPart = (_r[i+1]*_r[i+1]*_r[i+1] - _r[i]*_r[i]*_r[i])/3;
    double muPart = _mu[j+1] - _mu[j];
    double phiPart = _phi[k+1] - _phi[k];
    return rPart*muPart*phiPart;
}