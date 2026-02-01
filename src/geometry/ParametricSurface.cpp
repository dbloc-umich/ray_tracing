#include "ParametricSurface.h"
#include "GaussKronrod.h"
#include "NewtonSolver.h"
#include "ProjectedNewton.h"

ParametricSurface::ParametricSurface(property_ptr prop, const Eigen::Matrix3d& M,
                                     const Eigen::Vector3d& dr, std::shared_ptr<Material> mat):
    FunctionShape(prop, M, dr, mat)
{
    if (!_only_scaling){
        _lower = Eigen::VectorXd(3);
        _upper = Eigen::VectorXd(3);
        computeExtrema();
    }
    if (_sigma[0]/_sigma[2] - 1 > Shape::eps) computeSurfaceArea();
}

bool ParametricSurface::surfaceContains(const Eigen::Vector3d& p) const noexcept{
    // Least square problem to solve for
    auto f = [this, &p](const Eigen::Vector2d& u) -> Eigen::Vector3d { return p - r(u); };
    auto Jf = [this](const Eigen::Vector2d& u){
        Eigen::Matrix<double, 3, 2> jac;
        jac.col(0) = ru(u);
        jac.col(1) = rv(u);
        return jac;
    };
    Eigen::Vector2d u{(_prop->_u0+_prop->_u1)/2, (_prop->_v0+_prop->_v1)/2};
    NewtonSolver<2, 3> solver(f, Jf);
    auto status = solver.solve(u);
    switch (status){
        case (NLStatus::Success):
            return (_prop->_u0 <= u[0] && _prop->_u1 >= u[0] && _prop->_v0 <= u[1] && _prop->_v1 >= u[1]
                    && f(u).squaredNorm() <= Shape::eps * Shape::eps);
        default:
            return false;
    }
}

bool ParametricSurface::encloses(const Eigen::Vector3d& p) const noexcept{
    if (!_prop->_isClosedSurface) return false; // Checking if a point is inside is only defined for closed surfaces
    
    UnitVector3d dir(p - r({_prop->_u0, _prop->_v0}));
    return distanceToSurface(p, dir) > Shape::eps && distanceToSurface(p, -dir) > Shape::eps;
}

bool ParametricSurface::overlaps(const Shape& other) const noexcept{
    return false;
}

double ParametricSurface::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{
    auto f = [this, &p, &dir](const Eigen::Vector3d& s) -> Eigen::Vector3d
        { return p + dir.value()*s[2] - r({s[0], s[1]}); };
    auto Jf = [this, &dir](const Eigen::Vector3d& s){
        Eigen::Matrix3d jac;
        jac.col(0) = -ru({s[0], s[1]});
        jac.col(1) = -rv({s[0], s[1]});
        jac.col(2) = dir.value();
        return jac;
    };
    Eigen::Vector3d s{(_prop->_u0+_prop->_u1)/2, (_prop->_v0+_prop->_v1)/2, xMax()-xMin()};
    NewtonSolver<3> solver(f, Jf);
    auto status = solver.solve(s);
    switch (status){
        case (NLStatus::Success):
            if (_prop->_u0 < s[0] && _prop->_u1 > s[0] && _prop->_v0 < s[1] && _prop->_v1 > s[1])
                return std::numeric_limits<double>::quiet_NaN();
            if (s[2] > Shape::eps) return s[2];
            if (std::abs(s[2]) <= Shape::eps) return 0.0;
            return std::numeric_limits<double>::quiet_NaN();
        default:
            return std::numeric_limits<double>::quiet_NaN();
    }
}

UnitVector3d ParametricSurface::normal(const Eigen::Vector3d& pos) const{
    return UnitVector3d(nullptr);
}

void ParametricSurface::computeExtrema(){
    for (std::size_t i = 0; i < 3; i++){
        // Check the corners first      
        auto x = [this, i](const Eigen::Vector2d& u){ return _M.row(i).dot(_prop->r(u)); }; // Spatial variable to optimize    
        double xMin = std::numeric_limits<double>::max();
        double xMax = std::numeric_limits<double>::lowest();
        for (double u: {_prop->_u0, _prop->_u1}){
            for (double v: {_prop->_v0, _prop->_v1}){
                double val = x({u, v});
                if (val < xMin) xMin = val;
                if (val > xMax) xMax = val;
            }
        }

        auto grad = [this, i](const Eigen::Vector2d& u) -> Eigen::Vector2d
        {
            Eigen::Matrix<double, 3, 2> ruv;
            ruv.col(0) = _prop->ru(u);
            ruv.col(1) = _prop->rv(u);
            return (_M.row(i) * ruv).transpose();
        };
        auto H = [this, i](const Eigen::Vector2d& u) -> Eigen::Matrix2d
            { return _M(i,0)*_prop->_Hx(u) + _M(i,1)*_prop->_Hy(u) + _M(i,2)*_prop->_Hz(u); };

        // Sample values of [u, v] to find all critical points
        std::vector<Eigen::Vector2d> critPts;
        constexpr std::size_t Nu = 4;
        constexpr std::size_t Nv = 4;
        // Using Newton's method to find all possible critical points
        for (std::size_t i = 0; i < Nu; i++){
            for (std::size_t j = 0; j < Nv; j++){
                // Minimize
                Eigen::Vector2d argmin{_prop->_u0+(_prop->_u1-_prop->_u0)/(Nu-1)*i, _prop->_v0+(_prop->_v1-_prop->_v0)/(Nv-1)*j};
                Eigen::Vector2d argmax(argmin);
                ProjectedNewton<2> optimizer(x, {_prop->_u0, _prop->_v0}, {_prop->_u1, _prop->_v1}, grad, H);
                auto status = optimizer.minimize(argmin);
                if (status == COStatus::Success) xMin = std::min(xMin, x(argmin));

                // Maximize
                optimizer.setFunction([&x](const auto& u){ return -x(u); },
                                      [&grad](const auto& u) -> Eigen::Vector2d { return -grad(u); },
                                      [&H](const auto& u) -> Eigen::Matrix2d { return -H(u); });
                status = optimizer.minimize(argmax);
                if (status == COStatus::Success) xMax = std::max(xMax, x(argmax));
            }
        }

        _lower[i] = xMin + _dr[i];
        _upper[i] = xMax + _dr[i];
    }
}

void ParametricSurface::computeSurfaceArea(){
    auto dS = [this](const Eigen::Vector2d& u){
        Eigen::Vector3d v = _prop->ru(u).cross(_prop->rv(u));
        v[0] *= _sigma[1]*_sigma[2];
        v[1] *= _sigma[0]*_sigma[2];
        v[2] *= _sigma[0]*_sigma[1];
        return v.norm();
    };
    
    // Gauss-Kronrod integration with 15 points
    double uub = _prop->_u0 + (_prop->_u1-_prop->_u0)/_prop->_uSym; // upper bound in u
    double vub = _prop->_v0 + (_prop->_v1-_prop->_v0)/_prop->_vSym; // upper bound in v
    GaussKronrod<2> GL(7); // Gauss-Legendre integration with 5 points
    GaussKronrod<2>::IntegrationDomain D{_prop->_u0, uub, _prop->_v0, vub};
    _surfaceArea = GL.integrate(dS, D) * (_prop->_uSym*_prop->_vSym);
}

std::ostream& ParametricSurface::print(std::ostream& os) const noexcept{
    os << "Parametric surface defined by function stored at " << _prop.get();
    os << ", transformed by \n" << _M;
    os << ", and translated by [" << _dr[0] << ", " << _dr[1] << ", " << _dr[2] << "].";
    return os; 
}