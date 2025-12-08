#include "Ellipsoid.h"
#include "Box.h"
#include "Constants.h"
#include "NonlinearSolvers.h"
#include "Sphere.h"

#include "Eigen/Eigenvalues"

#include <cassert>
#include <cmath>


Ellipsoid::Ellipsoid(const Eigen::Vector3d& pt, double a, double b, double c, std::shared_ptr<Material> mat):
    Shape(mat),
    _center(pt),
    _M(), // empty; diagonal matrix with entries 1/a^2, 1/b^2, and 1/c^2
    _axes(), // empty; primary axes are i-hat, j-hat, and k-hat 
    _length(a, b, c)
{
    if (a <= 0 || b <= 0 || c <= 0) throw std::invalid_argument("ERROR: All semi-axes must be positive.");
}

Ellipsoid::Ellipsoid(const Eigen::Vector3d& pt, const Eigen::Matrix3d& M, std::shared_ptr<Material> mat):
    Shape(mat),
    _center(pt),
    _M(M)
{
    if (!M.isApprox(M.transpose())) throw std::invalid_argument("ERROR: The input matrix must be symmetric positive definite.");
    Eigen::SelfAdjointEigenSolver<std::decay_t<decltype(M)>> es(M);
    if (es.info() == Eigen::NoConvergence) throw std::runtime_error("ERROR: Eigen decomposition failed.");
    
    _axes = es.eigenvectors();
    _length = es.eigenvalues();
    if ((_length <= 0).any()) throw std::invalid_argument("ERROR: The input matrix must be symmetric positive definite.");
    _length = 1.0/(_length.sqrt());
}

double Ellipsoid::semiAxis(Eigen::Index i) const noexcept{
    assert(i < 3);
    return _length[i];
}

UnitVector3d Ellipsoid::principalAxis(Eigen::Index i) const noexcept{
    assert(i < 3);
    if (_axes.size() == 0){
        if (i == 0) return UnitVector3d(1, 0, 0);
        if (i == 1) return UnitVector3d(0, 1, 0);
        if (i == 2) return UnitVector3d(0, 0, 1);
    }
    Eigen::Vector3d axis = _axes.col(i).eval();
    return axis;
}

void Ellipsoid::setSemiAxis(Eigen::Index i, double L){
    assert(i < 3);
    if (L <= 0.0) throw std::invalid_argument("ERROR: Any semi-axis must be positive.");
    _length[i] = L;
    if (_M.size() != 0){
        // Recalculates the matrix
        Eigen::DiagonalMatrix<double, 3> Lambda{1.0 / _length / _length};
        _M = _axes * Lambda * _axes.transpose();
    }
}

double Ellipsoid::xMin() const noexcept{
    if (_M.size() == 0) return _center[0] - _length[0];
    Eigen::Vector3d v = _length * _axes.row(0).transpose().array();
    return _center[0] - std::sqrt(v.dot(v));
}

double Ellipsoid::xMax() const noexcept{
    if (_M.size() == 0) return _center[0] + _length[0];
    Eigen::Vector3d v = _length * _axes.row(0).transpose().array();
    return _center[0] + std::sqrt(v.dot(v));
}

double Ellipsoid::yMin() const noexcept{
    if (_M.size() == 0) return _center[1] - _length[1];
    Eigen::Vector3d v = _length * _axes.row(1).transpose().array();
    return _center[1] - std::sqrt(v.dot(v));
}

double Ellipsoid::yMax() const noexcept{
    if (_M.size() == 0) return _center[1] + _length[1];
    Eigen::Vector3d v = _length * _axes.row(1).transpose().array();
    return _center[1] + std::sqrt(v.dot(v));
}

double Ellipsoid::zMin() const noexcept{
    if (_M.size() == 0) return _center[2] - _length[2];
    Eigen::Vector3d v = _length * _axes.row(2).transpose().array();
    return _center[2] - std::sqrt(v.dot(v));
}

double Ellipsoid::zMax() const noexcept{
    if (_M.size() == 0) return _center[2] + _length[2];
    Eigen::Vector3d v = _length * _axes.row(2).transpose().array();
    return _center[2] + std::sqrt(v.dot(v));
}

double Ellipsoid::surfaceArea() const noexcept{
    auto norm_normal = [this](double u, double v){
        auto cos2v = std::cos(v);
        cos2v *= cos2v;
        return std::sqrt(xdotMy(Eigen::Vector3d{1.0-u*u, 1.0-u*u, u*u}, Eigen::Vector3d{cos2v, 1-cos2v, 1}));
    };

    // Gauss-Legendre integration with 3 points
    double S = 0.0;
    const Eigen::Array3d roots{-std::sqrt(3.0/5), 0.0, std::sqrt(3.0/5)};
    const Eigen::Array3d weights{5.0/9, 8.0/9, 5.0/9};
    for (Eigen::Index i = 0; i < 3; i++){
        // Integrating over u in [0, 1]
        double u = 0.5 * (roots[i] + 1);
        double wu = 0.5 * weights[i];
        for (Eigen::Index j = 0; j < 3; j++){
            // Integrating over v in [0, PI/2]
            double v = mconst::pi/4 * (roots[j] + 1);
            double wv = mconst::pi/4 * weights[j];
            S += wu * wv * norm_normal(u, v);
        }
    }
    return 8*_length[0]*_length[1]*_length[2]*S;
}
double Ellipsoid::volume() const noexcept{ return 4.0/3*mconst::pi*_length[0]*_length[1]*_length[2]; }

bool Ellipsoid::surfaceContains(const Eigen::Vector3d& p) const noexcept{
    Eigen::Vector3d dr = _center - p;
    if (_M.size() == 0) return std::abs(xdotMy(dr, dr) - 1) <= Shape::eps;
    return std::abs(dr.dot(_M*dr) - 1) <= Shape::eps;
}

bool Ellipsoid::encloses(const Eigen::Vector3d& p) const noexcept{
    Eigen::Vector3d dr = _center - p;
    if (_M.size() == 0) return xdotMy(dr, dr) - 1 < -Shape::eps;
    return dr.dot(_M*dr) - 1 < -Shape::eps;
}

bool Ellipsoid::overlaps(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object
    
    const Sphere* sphere = dynamic_cast<const Sphere*>(&other);
    const Ellipsoid* ellipsoid = dynamic_cast<const Ellipsoid*>(&other);
    if (sphere || ellipsoid){
        Eigen::Matrix3d Binv;
        Eigen::Vector3d dr = _center;
        if (sphere){
            Binv = Eigen::Matrix3d::Identity() * sphere->radius() * sphere->radius();
            dr -= sphere->center();
        } else{
            Binv = Eigen::DiagonalMatrix<double, 3>(ellipsoid->_length[0] * ellipsoid->_length[0],
                                                    ellipsoid->_length[1] * ellipsoid->_length[1],
                                                    ellipsoid->_length[2] * ellipsoid->_length[2]);
            dr -= ellipsoid->center();
            if (ellipsoid->_M.size() != 0) Binv = ellipsoid->_axes * Binv * ellipsoid->_axes.transpose();
        }
        if (dr.squaredNorm() <= Shape::eps) return true; // centers overlap

        auto K = [this, &dr, &Binv](double lambda){
            Eigen::Matrix3d M = Binv/(1.0-lambda);
            Eigen::Matrix3d Ainv = Eigen::DiagonalMatrix<double, 3>(_length[0]*_length[0], _length[1]*_length[1], _length[2]*_length[2]);
            if (_M.size() != 0) Ainv = _axes * Ainv * _axes.transpose();
            M += Ainv/lambda;
            Eigen::PartialPivLU<Eigen::Ref<Eigen::Matrix3d>> lu(M);
            return 1.0 - dr.dot(lu.solve(dr));
        };
        try{
            double lambda = newton(K, 0.5); // Find a root of the K(lambda)
            if (lambda > 0.0 && lambda < 1.0) return false; // root is found on (0, 1), no overlapping
            return true; // no root found on (0, 1)
        } catch(const NewtonSingularityError& ex){
            return false; // double root is found
        } catch(const NewtonMaxIterationsError& ex){
            return true; // no root is found
        }
    }

    if (const Box* box = dynamic_cast<const Box*>(&other)){
        if (box->encloses(_center)) return true; // Ellipsoid center check
        auto p0 = box->lowerVertex();
        auto p1 = box->upperVertex();
        if (this->encloses((p0+p1)/2)) return true; // Box center check

        // Checks closest vertices
        std::vector<Eigen::Vector3d> vertices;
        for (auto x: {p0[0], p1[0]}){
            for (auto y: {p0[1], p1[1]}){
                for (auto z: {p0[2], p1[2]}){
                    vertices.emplace_back(x, y, z);
                }
            }
        }

        auto removalCriteria = [this, &p0, &p1](const auto& v){
            for (Eigen::Index i = 0; i < 3; i++){
                if (_center[i] < p0[i] && v[i] == p0[i]) return false;
                if (_center[i] > p1[i] && v[i] == p1[i]) return false;
            }
            return true;
        };
        vertices.erase(std::remove_if(vertices.begin(), vertices.end(), removalCriteria), vertices.end());

        for (auto vertex: vertices){
            if (this->encloses(vertex)) return true;
        }
        
        // Checks closest edges
        for (auto it1 = vertices.begin(); it1 != vertices.end(); it1++){
            for (auto it2 = it1+1; it2 != vertices.end(); it2++){
                Eigen::Vector3d diff = *it2 - *it1;
                // Check if the two points form an edge, then intersection check
                if ((diff.array() != 0).count() == 1 && this->intersects(*it1, *it2)) return true;
            }
        }

        // If all four vertices of a face is checked, remove a fourth vertex to avoid future redundant face checks
        Eigen::Array3d faces;
        faces.fill(std::numeric_limits<double>::quiet_NaN());
        for (auto it1 = vertices.begin(); it1 != vertices.end();){
            Eigen::Array<std::size_t, 3, 1> counts{0, 0, 0}; // the number of other vertices that share its coordinate
            for (auto it2 = it1+1; it2 != vertices.end(); ++it2){
                for (Eigen::Index i = 0; i < 3; i++){
                    if ((*it1)[i] == (*it2)[i]) counts[i]++;
                }
            }

            bool toBeRemoved = false;
            for (Eigen::Index i = 0; i < 3; i++){
                if (counts[i] == 3 && faces[i] != (*it1)[i]){
                    faces[i] = (*it1)[i];
                    toBeRemoved = true;
                }
            }
            if (toBeRemoved) vertices.erase(it1);
            else ++it1;
        }

        // Checks closest faces
        for (auto it1 = vertices.begin(); it1 != vertices.end(); ++it1){
            for (auto it2 = it1+1; it2 != vertices.end(); ++it2){
                for (auto it3 = it2+1; it3 != vertices.end(); ++it3){
                    for (Eigen::Index i = 0; i < 3; i++){
                        // Confirm that the three points are on the same face
                        if ((*it1)[i] == (*it2)[i] && (*it1)[i] == (*it3)[i]){                    
                            // Find the point at the right angle, and then intersection check
                            Eigen::Vector3d d21 = *it2 - *it1;
                            Eigen::Vector3d d31 = *it3 - *it1;
                            Eigen::Vector3d d32 = *it3 - *it2;
                            if (d21.dot(d31) == 0.0){
                                // *it1 is at the right angle
                                if (this->intersects(*it1, *it2, *it3)) return true;
                            }
                            if (d21.dot(d32) == 0.0){
                                // *it2 is at the right angle
                                if (this->intersects(*it2, *it3, *it1)) return true;
                            }
                            if (this->intersects(*it3, *it1, *it2)) return true; // *it3 is at the right angle
                            break;
                        }
                    }
                }
            }
        }

        return false;
    }
    return other.overlaps(*this);
}

double Ellipsoid::distanceToSurface(const Eigen::Vector3d& p, const UnitVector3d& dir) const noexcept{
    Eigen::Vector3d dr = p - _center;
    Eigen::Vector3d Omega = dir.value();
    double A, B, C;
    if (_M.size() == 0){
        A = xdotMy(Omega, Omega);
        B = xdotMy(Omega, dr);
        C = xdotMy(dr, dr) - 1;
    } else{
        A = Omega.dot(_M*Omega);
        B = Omega.dot(_M*dr);
        C = dr.dot(_M*dr) - 1;
    }

    if (std::abs(C) < Shape::eps){ // on the surface
        if (B >= 0.0) return 0.0; // leaves the surface
        return -2*B/A; // enters the surface and leaves at another location
    }

    double discr = B*B - A*C;
    if (discr < 0) return NAN; // particle is outside and will never enters the Ellipsoid
    double sp, sm; // two roots of the solutions
    if (B < 0){
        sp = (-B + std::sqrt(discr))/A;
        sm = C/A/sp; 
    } else{
        sm = (-B - std::sqrt(discr))/A;
        sp = C/A/sm;
    }

    if (C > 0) return (sp > 0 ? std::min(sp,sm) : NAN); // particle is outside, the roots have the same sign
    return (sp > 0 ? sp : sm); // the roots have opposite sign, take the positive one - particle is inside
}

UnitVector3d Ellipsoid::normal(const Eigen::Vector3d& pos) const{
    if (!surfaceContains(pos)) throw std::invalid_argument("ERROR: Eigen::Vector3d is not on the surface.");
    Eigen::Vector3d dr = pos - _center;
    if (_M.size() == 0){
        dr[0] /= (_length[0]*_length[0]);
        dr[1] /= (_length[1]*_length[1]);
        dr[2] /= (_length[2]*_length[2]);
    } else dr = _M*dr;
    return dr;
}

std::ostream& Ellipsoid::print(std::ostream& os) const noexcept{
    if (_M.size() == 0){
        if (_center[0] == 0.0) os << "x^2";
        else{
            os << "(x ";
            if (_center[0] > 0) os << "- " << _center[0] << ")^2";
            else os << "+ " << -_center[0] << ")^2";
        }
        os << "/" << _length[0]*_length[0] << " + ";

        if (_center[1] == 0.0) os << "y^2";
        else{
            os << "(y ";
            if (_center[1] > 0) os << "- " << _center[1] << ")^2";
            else os << "+ " << -_center[1] << ")^2";
        }
        os << "/" << _length[1]*_length[1] << " + ";

        if (_center[2] == 0.0) os << "z^2";
        else{
            os << "(z ";
            if (_center[2] > 0) os << "- " << _center[2] << ")^2";
            else os << "+ " << -_center[2] << ")^2";
        }
        os << "/" << _length[2]*_length[2] << " = 1";
    } else{
        os << "(r - [" << _center[0] << ", " << _center[1] << ", " << _center[2] << "]^T)^T";
        os << _M;
        os << "(r - [" << _center[0] << ", " << _center[1] << ", " << _center[2] << "]^T) = 1";
    }
    return os;
}

double Ellipsoid::xdotMy(const Eigen::Vector3d& x, const Eigen::Vector3d& y) const noexcept{
    // Compute the dot product x^T * My for an axis-aligned ellipsoid
    return x[0]*y[0]/(_length[0]*_length[0]) + x[1]*y[1]/(_length[1]*_length[1]) + x[2]*y[2]/(_length[2]*_length[2]);
}

bool Ellipsoid::intersects(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1) const noexcept{
    // Make sure p0 and p1 are both outside of the ellipsoid
    Eigen::Vector3d p1p0 = p1 - p0;
    Eigen::Vector3d p0r0 = p0 - _center;
    double A, B, C;
    if (_M.size() == 0){
        A = xdotMy(p1p0, p1p0);
        B = xdotMy(p1p0, p0r0);
        C = xdotMy(p0r0, p0r0);
    }
    else{
        A = p1p0.dot(_M*p1p0);
        B = p1p0.dot(_M*p0r0);
        C = p0r0.dot(_M*p0r0);
    }
    double t = -B/A;
    return t >= -1 && t <= 0 && A*t*t + 2*B*t + C < 1;
}

bool Ellipsoid::intersects(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) const noexcept{
    // Make sure p0, p1, p2 are all outside of the ellipsoid, and the line p0-p1 is perpendicular to p0-p2.
    Eigen::Vector3d p1p0 = p1 - p0;
    Eigen::Vector3d p2p0 = p2 - p0;
    Eigen::Vector3d p0r0 = p0 - _center;
    double A, B, C, D, E, F;
    if (_M.size() == 0){
        A = xdotMy(p1p0, p1p0);
        B = xdotMy(p1p0, p2p0);
        C = xdotMy(p2p0, p2p0);
        D = xdotMy(p1p0, p0r0);
        E = xdotMy(p2p0, p0r0);
        F = xdotMy(p0r0, p0r0);
    }
    else{
        A = p1p0.dot(_M*p1p0);
        B = p1p0.dot(_M*p2p0);
        C = p2p0.dot(_M*p2p0);
        D = p1p0.dot(_M*p0r0);
        E = p2p0.dot(_M*p0r0);
        F = p0r0.dot(_M*p0r0);
    }
    double BBAC = B*B - A*C;
    if (std::abs(BBAC) <= Shape::eps) return false;
    double u = (C*D - B*E)/BBAC;
    double v = (A*E - B*D)/BBAC;
    return u >= 0 && u <= 1 && v >= 0 && v <= 1 && A*u*u + 2*B*u*v + C*v*v + 2*D*u + 2*E*v + F < 1;
}