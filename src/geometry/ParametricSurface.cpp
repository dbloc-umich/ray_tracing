#include "ParametricSurface.h"
#include "NonlinearSolver.h"
#include "Optimizer.h"

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
    auto f = [this, &p](const Eigen::Vector2d& u){ return p - r(u); };
    auto Jf = [this](const Eigen::Vector2d& u){
        Eigen::Matrix<double, 3, 2> jac;
        jac.col(0) = ru(u);
        jac.col(1) = rv(u);
        return jac;
    };
    Eigen::Vector2d u{(_prop->_u0+_prop->_u1)/2, (_prop->_v0+_prop->_v1)/2};
    try{
        newton(f, Jf, u);
        if (_prop->_u0 <= u[0] && _prop->_u1 >= u[0] && _prop->_v0 <= u[1] && _prop->_v1 >= u[1]
            && f(u).squaredNorm() <= Shape::eps) return true;
        return false;
    } catch (const std::runtime_error&){
        return false;
    }
}

bool ParametricSurface::encloses(const Eigen::Vector3d& p) const noexcept{
    if (!_prop->_isClosedSurface) return false;
        // Checking if a point is inside is only defined for closed surfaces
    
    UnitVector3d dir(p - r({_prop->_u0, _prop->_v0}));
    return distanceToSurface(p, dir) > Shape::eps && distanceToSurface(p, -dir) > Shape::eps;
}

bool ParametricSurface::overlaps(const Shape& other) const noexcept{
    return false;
}

double ParametricSurface::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{
    auto f = [this, &p, &dir](const Eigen::Vector3d& s){ return p + dir.value()*s[2] - r({s[0], s[1]}); };
    auto Jf = [this, &dir](const Eigen::Vector3d& s){
        Eigen::Matrix3d jac;
        jac.col(0) = -ru({s[0], s[1]});
        jac.col(1) = -rv({s[0], s[1]});
        jac.col(2) = dir.value();
        return jac;
    };
    Eigen::Vector3d s{(_prop->_u0+_prop->_u1)/2, (_prop->_v0+_prop->_v1)/2, xMax()-xMin()};
    try{
        newton(f, Jf, s);
        if (_prop->_u0 <= s[0] && _prop->_u1 >= s[0] && _prop->_v0 <= s[1] && _prop->_v1 >= s[1])
            return std::numeric_limits<double>::quiet_NaN();
        if (s[2] > Shape::eps) return s[2];
        if (std::abs(s[2]) <= Shape::eps) return 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    } catch (const std::runtime_error&){
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
                Eigen::Vector2d argmin{_prop->_u0+(_prop->_u1-_prop->_u0)/(Nu-1)*i, _prop->_v0+(_prop->_v1-_prop->_v0)/(Nv-1)*j};
                Eigen::Vector2d argmax = argmin;
                try{
                    projectedNewton(x, argmin, {_prop->_u0, _prop->_v0}, {_prop->_u1, _prop->_v1}, grad, H);
                    double val = x(argmin);
                    if (val < xMin) xMin = val;

                    projectedNewton([&x](const auto& u){ return -x(u); },
                                    argmax, {_prop->_u0, _prop->_v0}, {_prop->_u1, _prop->_v1},
                                    [&grad](const auto& u) -> Eigen::Vector2d { return -grad(u); },
                                    [&H](const auto& u) -> Eigen::Matrix2d { return -H(u); });
                    val = x(argmax);
                    if (val > xMax) xMax = val;
                } catch(const std::runtime_error& ex){}
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
    
    // Gauss-Legendre integration with 3 points
    _surfaceArea = 0.0;
    double uub = _prop->_u0 + (_prop->_u1-_prop->_u0)/_prop->_uSym; // upper bound in u
    double vub = _prop->_v0 + (_prop->_v1-_prop->_v0)/_prop->_vSym; // upper bound in v
    const Eigen::Array3d roots{-std::sqrt(3.0/5), 0.0, std::sqrt(3.0/5)};
    const Eigen::Array3d weights{5.0/9, 8.0/9, 5.0/9};
    
    for (Eigen::Index i = 0; i < 3; i++){
        double u = ((uub-_prop->_u0)*roots[i] + uub + _prop->_u0)/2;
        double wu = (uub-_prop->_u0)/2 * weights[i];
        for (Eigen::Index j = 0; j < 3; j++){
            double v = ((vub-_prop->_v0)*roots[j] + vub + _prop->_v0)/2;
            double wv = (vub-_prop->_v0)/2 * weights[j];
            _surfaceArea += wu * wv * dS({u,v});
        }
    }
    _surfaceArea *= (_prop->_uSym*_prop->_vSym);
}

std::ostream& ParametricSurface::print(std::ostream& os) const noexcept{
    os << "Parametric surface defined by function stored at " << _prop.get();
    os << ", transformed by \n" << _M;
    os << ", and translated by [" << _dr[0] << ", " << _dr[1] << ", " << _dr[2] << "].";
    return os; 
}