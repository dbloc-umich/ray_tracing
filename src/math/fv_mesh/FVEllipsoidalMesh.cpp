#include "FVEllipsoidalMesh.h"
#include "Constants.h"
#include "Ellipsoid.h"
#include "GaussKronrod.h"
#include "NewtonSolver.h"
#include "Eigen/Cholesky"

FVEllipsoidalMesh::FVEllipsoidalMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, double a, double b, double c):
    FVSpatialMesh(x, y, z),
    _a(std::max({a, b, c})),
    _b(a+b+c - std::max({a, b, c}) - std::min({a, b, c})),
    _c(std::min({a, b, c})),
    _r(_axes[0]),
    _mu(_axes[1]),
    _phi(_axes[2])
{
    if (_r[0] < 0 || _r[_r.size()-1] > 1) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_mu[0] < -1 || _mu[_mu.size()-1] > 1) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_phi[0] < -mconst::pi || _phi[_phi.size()-1] > mconst::pi) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_a <= 0 || _b <= 0 || _c <= 0) throw std::invalid_argument("ERROR: Non-positive axis length.");
}

FVEllipsoidalMesh::FVEllipsoidalMesh(const std::array<Eigen::ArrayXd, 3>& axes, double a, double b, double c):
    FVSpatialMesh(axes),
    _a(std::max({a, b, c})),
    _b(a+b+c - std::max({a, b, c}) - std::min({a, b, c})),
    _c(std::min({a, b, c})),
    _r(_axes[0]),
    _mu(_axes[1]),
    _phi(_axes[2])
{
    if (_r[0] < 0 || _r[_r.size()-1] > 1) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_mu[0] < -1 || _mu[_mu.size()-1] > 1) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_phi[0] < -mconst::pi || _phi[_phi.size()-1] > mconst::pi) throw std::invalid_argument("ERROR: Invalid boundary values.");
    if (_a <= 0 || _b <= 0 || _c <= 0) throw std::invalid_argument("ERROR: Non-positive axis length.");
}

FVEllipsoidalMesh::FVEllipsoidalMesh(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, const Ellipsoid& ell):
    FVEllipsoidalMesh(x, y, z, ell.semiAxis(0), ell.semiAxis(1), ell.semiAxis(2))
{}

FVEllipsoidalMesh::FVEllipsoidalMesh(const std::array<Eigen::ArrayXd, 3>& axes, const Ellipsoid& ell):
    FVEllipsoidalMesh(axes, ell.semiAxis(0), ell.semiAxis(1), ell.semiAxis(2))
{}

Eigen::Vector3d FVEllipsoidalMesh::cartesian(double r, double mu, double phi) const noexcept{
    return {_a*r*mu, _b*r*std::sqrt(1-mu*mu)*std::cos(phi), _c*r*std::sqrt(1-mu*mu)*std::sin(phi)};
}

Eigen::Vector3d FVEllipsoidalMesh::coordinates(double x, double y, double z) const noexcept{
    x /= _a;
    y /= _b;
    z /= _c;
    double r = std::sqrt(x*x + y*y + z*z);
    return {r, x/r, std::atan2(z,y)};
}

Eigen::Vector3d FVEllipsoidalMesh::basisVector(Eigen::Index varID, double r, double mu, double phi) const noexcept{
    // varID = 0, 1, 2 referring to which basis vector to return
    if (varID == 0) return {_a*mu, _b*std::sqrt(1-mu*mu)*std::cos(phi), _c*std::sqrt(1-mu*mu)*std::sin(phi)};
    if (varID == 1) return {_a*r, -_b*r*mu/std::sqrt(1-mu*mu)*std::cos(phi), -_c*r*mu/std::sqrt(1-mu*mu)*std::sin(phi)};
    return {0, -_b*r*std::sqrt(1-mu*mu)*std::sin(phi), _c*r*std::sqrt(1-mu*mu)*std::cos(phi)};
}

UnitVector3d FVEllipsoidalMesh::normalVector(Eigen::Index varID, double r, double mu, double phi) const noexcept{
    if (varID == 0){
        if (r == 0) return UnitVector3d(nullptr);
        return basisVector(2, r, mu, phi).cross(basisVector(1, r, mu, phi));
    }
    if (varID == 1){
        if (mu == 1 || mu == -1) return UnitVector3d(nullptr);
        return basisVector(0, r, mu, phi).cross(basisVector(2, r, mu, phi));
    }
    return basisVector(1, r, mu, phi).cross(basisVector(0, r, mu, phi));
}

double FVEllipsoidalMesh::gradientDotN(double dudq, Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const{
    double r   = surfID == 0 ? _r[i]   : (surfID == 1 ? _r[i+1]   : (_r[i+1]+_r[i])/2);
    double mu  = surfID == 2 ? _mu[j]  : (surfID == 3 ? _mu[j+1]  : (_mu[j+1]+_mu[j])/2);
    double phi = surfID == 4 ? _phi[k] : (surfID == 5 ? _phi[k+1] : (_phi[k+1]+_phi[k])/2);
    
    Eigen::Vector3d er   = basisVector(0, r, mu, phi);
    Eigen::Vector3d emu  = basisVector(1, r, mu, phi);
    Eigen::Vector3d ephi = basisVector(2, r, mu, phi);

    Eigen::Matrix3d J; // Jacobian matrix
    J.col(0) = er;
    J.col(1) = emu;
    J.col(2) = ephi;

    Eigen::Matrix3d G; // Gram matrix
    G(0,0) = er.dot(er);
    G(0,1) = er.dot(emu);
    G(0,2) = er.dot(ephi);
    G(1,0) = G(0,1);
    G(1,1) = emu.dot(emu);
    G(1,2) = emu.dot(ephi);
    G(2,0) = G(0,2);
    G(2,1) = G(1,2);
    G(2,2) = ephi.dot(ephi);
    Eigen::LDLT<Eigen::Ref<Eigen::Matrix3d>> ldlt(G); // LDL^T decomposition

    // Function to find the gradient
    auto gradFunc = [&, this](const Eigen::Vector3d& g){
        Eigen::Vector3d rhs;
        rhs[0] = (surfID == 0 || surfID == 1) ? dudq : g.dot(cartesian(_r[i+1], mu, phi) - cartesian(_r[i], mu, phi))/(_r[i+1]-_r[i]);
        rhs[1] = (surfID == 2 || surfID == 3) ? dudq : g.dot(cartesian(r, _mu[j+1], phi) - cartesian(r, _mu[j], phi))/(_mu[j+1]-_mu[j]);
        rhs[2] = (surfID == 4 || surfID == 5) ? dudq : g.dot(cartesian(r, mu, _phi[k+1]) - cartesian(r, mu, _phi[k]))/(_phi[k+1]-_phi[k]);
        return g - J*ldlt.solve(rhs);
    };
    NewtonSolver<3> newton(gradFunc); // should be a linear system of equations, but use Newton's regardless
    newton.setFTol(1e-9);
    Eigen::Vector3d grad{dudq, dudq, dudq};
    auto status = newton.solve(grad);
    if (status == NLStatus::Success){
        UnitVector3d n = normalVector(surfID/2, r, mu, phi);
        return grad.dot(n.value());
    }
    throw std::runtime_error("ERROR: Unable to solve for the gradient.");
}

double FVEllipsoidalMesh::area(Eigen::Index i, Eigen::Index j, Eigen::Index k, Eigen::Index surfID) const noexcept{
    double del = _c/_b;
    double eps = _c/_a;
    double S;
    if (surfID <= 1){
        // surfaces where r is constant
        double r = surfID == 0 ? _r[i] : _r[i+1];
        if (r == 0.0) return 0;
        double rPart = r*r;
        auto integrand = [this, k, del, eps](double mu){
            double K = (1-del*del)*(1-mu*mu) / (1 - (1-eps*eps)*mu*mu);
            K = std::sqrt(K);
            double phiPart = std::ellint_2(K, mconst::pi/2 - _phi[k]) - std::ellint_2(K, mconst::pi/2 - _phi[k+1]);
            return std::sqrt(1 - (1*eps*eps)*mu*mu) * phiPart;
        };
        double muPart = GaussKronrod<>(7).integrate(integrand, {_mu[j], _mu[j+1]});
        S = rPart * muPart;
    } else if (surfID <= 3){
        // surfaces where mu is constant
        double mu = surfID == 2 ? _mu[j] : _mu[j+1];
        if (mu == 1.0 || mu == -1.0) return 0;
        double rPart = (_r[i+1]*_r[i+1] - _r[i]*_r[i])/2;
        double muPart = std::sqrt((1-mu*mu) * (eps*eps + (1-eps*eps)*mu*mu));
        double K = mu*mu*(1-del*del) / (eps*eps + (1-eps*eps)*mu*mu);
        K = std::sqrt(K);
        double phiPart = std::ellint_2(K, mconst::pi/2 - _phi[k]) - std::ellint_2(K, mconst::pi/2 - _phi[k+1]);
        S = rPart * muPart * phiPart;
    } else{
        // surfaces where phi is constant
        double phi = surfID == 4 ? _phi[k] : _phi[k+1];
        double rPart = (_r[i+1]*_r[i+1] - _r[i]*_r[i])/2;
        double muPart = std::asin(_mu[j+1]) - std::asin(_mu[j]);
        double phiPart = std::sqrt(del*del + (1-del*del)*std::cos(phi)*std::cos(phi));
        S = rPart * muPart * phiPart;
    }
    return S*_a*_b;
}

double FVEllipsoidalMesh::volume(Eigen::Index i, Eigen::Index j, Eigen::Index k) const noexcept{
    double rPart = (_r[i+1]*_r[i+1]*_r[i+1] - _r[i]*_r[i]*_r[i])/3;
    double muPart = _mu[j+1] - _mu[j];
    double phiPart = _phi[k+1] - _phi[k];
    return _a*_b*_c*rPart*muPart*phiPart;
}

double FVEllipsoidalMesh::ghostLength(Eigen::Index varID) const noexcept{
    assert (varID == 1); // can only extend in the outer radial surface
    double r0 = _r[_r.size()-2]; // second to last boundary;
    return std::cbrt(2-r0); 
}