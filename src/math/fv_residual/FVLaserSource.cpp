#include "FVLaserSource.h"
#include "Constants.h"
#include "Ellipsoid.h"
#include "FVSpatialMesh.h"
#include "FVStateMesh.h"
#include "Ray.h"
#include "Sphere.h"

#include <random>
//#define MONITOR
#define MONITOR2

namespace{
    static constexpr double eps = 1.0e-9;
    static constexpr double IThreshold = 1.0e-9;
    static std::default_random_engine rng(64);
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
}

FVLaserSource::FVLaserSource(std::shared_ptr<Material> mat, Ray& ray, Shape* shape):
    FVResidual(mat, PropVariable::temperature),
    _ray(ray),
    _shape(shape),
    _I0(ray.intensity())
{
    assert(!_shape->surfaceContains(_ray.position()) && !_shape->encloses(_ray.position()));
}

Eigen::VectorXd FVLaserSource::computeResidual(const FVStateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nr   = mesh->axisSize(0)-1; // number of cells on the r axis
    Eigen::Index Nmu  = mesh->axisSize(1)-1; // number of cells on the mu axis
    Eigen::Index Nphi = mesh->axisSize(2)-1; // number of cells on the phi axis
    FVStateMesh q(mesh, 0);
    double a, b, c; // semi-axes of the ellipse
    bool isSpherical;
    if (Sphere* sph = dynamic_cast<Sphere*>(_shape)){
        a = sph->radius();
        b = sph->radius();
        c = sph->radius();
        isSpherical = true;
    } else if (Ellipsoid* ell = dynamic_cast<Ellipsoid*>(_shape)){
        a = std::max({ell->semiAxis(0), ell->semiAxis(1), ell->semiAxis(2)});
        c = std::min({ell->semiAxis(0), ell->semiAxis(1), ell->semiAxis(2)});
        b = ell->semiAxis(0) + ell->semiAxis(1) + ell->semiAxis(2) - a - c;
        isSpherical = false;
    } else throw std::invalid_argument("ERROR: Laser heat source has only been implemented for spherical and ellipsoidal meshes");

    Ray ray = _ray;
    double s = _shape->distanceToSurface(ray.position(), ray.direction());
    if (std::isnan(s)) return q.matrix(); // ray does not enter shape

    // Advances the ray forward to the medium surface
    Eigen::Vector3d omega = ray.direction().value();
    ray.setPoistion(ray.position() + omega*s);

    /*
    convert current ray position into the shape coordinates   
    */
    
    // Locates the first cell that the ray enters
    Eigen::Vector3d coord = mesh->coordinates(ray.position()[0], ray.position()[1], ray.position()[2]);
    Eigen::Index i = Nr-1; // only works for spherical and ellipsoidal coordinates
    Eigen::Index j = std::upper_bound(mesh->axis(1).cbegin(), mesh->axis(1).cend(), coord[1]) - mesh->axis(1).cbegin() - 1;
    Eigen::Index k = std::upper_bound(mesh->axis(2).cbegin(), mesh->axis(2).cend(), coord[2]) - mesh->axis(2).cbegin() - 1;

    // Refract the ray
    UnitVector3d normal = mesh->normalVector(0, coord[0], coord[1], ray.position()[2]);
    double n1 = 1.0;
    Material::PropVars vars{{_var, u(i,j,k)}};
    double n2 = _mat->computeProperty(Prop::refractiveIndex, vars);
    ray.setDirection(ray.direction().refract(normal, n1, n2));
    omega = ray.direction().value();

    while (true){
#ifdef MONITOR
        std::cout << "Ray is at " << ray.position().transpose() << ", which is " << mesh->coordinates(ray.position()).transpose() << " in spherical coordinates." << std::endl;
        std::cout << "Current direction is " << ray.direction() << std::endl;
        std::cout << "Current intensity is " << ray.intensity() << std::endl;
        std::cout << "Currently in cell (" << i << "," << j << "," << k << "): ";
        std::cout << "r in [" << mesh->axis(0)[i] << ", " << mesh->axis(0)[i+1] << "], ";
        std::cout << "mu in [" << mesh->axis(1)[j] << ", " << mesh->axis(1)[j+1] << "], ";
        std::cout << "phi in [" << mesh->axis(2)[k] << ", " << mesh->axis(2)[k+1] << "]" << std::endl;
#endif
#ifdef MONITOR2
        std::cout << "[" << ray.position().transpose() << "] on cell " << i << ", " << j << ", " << k << std::endl;
#endif
        s = std::numeric_limits<double>::max();
        std::size_t surfID = 6; // surfID of the next surface, 6 is an invalid ID

        // Find the exit surface and its distance
        for (std::size_t surf = 0; surf < 6; surf++){
            if (surf == 0 || surf == 1){
                // Check surfaces r = r_{i-1/2} and r_{i+1/2}
                if (surf == 0 && i == 0) continue; // can theoretically passes through the origin, but it's computationally impossible
                
                double r = surf == 0 ? mesh->axis(0)[i] : mesh->axis(0)[i+1];
                double A, B, C;
                if (isSpherical){
                    A = 1.0;
                    B = omega.dot(ray.position());
                    C = ray.position().squaredNorm() - r*r;
                } else{
                    Eigen::DiagonalMatrix<double, 3> M{1/(a*a), 1/(b*b), 1/(c*c)};
                    A = omega.dot(M*omega);
                    B = omega.dot(M*ray.position());
                    C = ray.position().dot(M*ray.position()) - r*r;
                }
                
                double discr = B*B - A*C;
                if (std::abs(discr) < eps) discr = 0;
                if (discr < 0) continue; // does not reach this surface
                double sp, sm; // two roots of the quadratic equation
                if (B < 0){
                    sp = (-B + std::sqrt(discr))/A;
                    sm = C/A/sp; 
                } else{
                    sm = (-B - std::sqrt(discr))/A;
                    sp = C/A/sm;
                }

                // Bound checking the roots
                for (auto dist: {sp, sm}){
                    Eigen::Vector3d dest = mesh->coordinates(ray.position() + omega*dist);
#ifdef MONITOR
                    std::cout << "On surface r = " << r << ", one possible distance is s = " << dist << ". ";
                    std::cout << "The destination coordinates are " << dest.transpose() << std::endl;
#endif
                    bool validS = dist > eps && dist < s;
                    bool validMu = dest[1] >= mesh->axis(1)[j] && dest[1] <= mesh->axis(1)[j+1];
                    bool validPhi = dest[2] >= mesh->axis(2)[k] && dest[2] <= mesh->axis(2)[k+1];
                    if (validS && validMu && validPhi){
                        s = dist;
                        surfID = surf;
#ifdef MONITOR
                        std::cout << "\tThis distance is accepted as the new best estimate." << std::endl;
#endif
                    }
#ifdef MONITOR
                    else std::cout << "\tThis distance is rejected." << std::endl;
#endif
                }
            } else if (surf == 2 || surf == 3){
                // Check surfaces mu = mu_{j-1/2} and mu_{j+1/2}
                if ((surf == 2 && j == 0) || (surf == 3 && j == Nmu-1)) continue; // can theoretically passes through the major axis, but it's computationally impossible
                
                double mu = surf == 2 ? mesh->axis(1)[j] : mesh->axis(1)[j+1];
                Eigen::DiagonalMatrix<double, 3> M{(1-mu*mu)/(a*a), -mu*mu/(b*b), -mu*mu/(c*c)};
                double A = omega.dot(M*omega);
                double B = omega.dot(M*ray.position());
                double C = ray.position().dot(M*ray.position());
                
                double discr = B*B - A*C;
                if (std::abs(discr) < eps) discr = 0;
                if (discr < 0) continue; // does not reach this surface
                double sp, sm; // two roots of the quadratic equation
                if (B < 0){
                    sp = (-B + std::sqrt(discr))/A;
                    sm = C/A/sp; 
                } else{
                    sm = (-B - std::sqrt(discr))/A;
                    sp = C/A/sm;
                }

                // Bound checking the roots
                for (auto dist: {sp, sm}){
                    Eigen::Vector3d dest = mesh->coordinates(ray.position() + omega*dist);
#ifdef MONITOR
                    std::cout << "On surface mu = " << mu << ", one possible distance is s = " << dist << ". ";
                    std::cout << "The destination coordinates are " << dest.transpose() << std::endl;
#endif
                    bool validS = dist > eps && dist < s;
                    bool validR = dest[0] >= mesh->axis(0)[i] && dest[0] <= mesh->axis(0)[i+1];
                    bool validMu = (ray.position()[0] + omega[0]*dist)*mu >= 0; // make sure that it maps to the surface mu and not -mu
                    bool validPhi = dest[2] >= mesh->axis(2)[k] && dest[2] <= mesh->axis(2)[k+1];
                    if (validS && validR && validMu && validPhi){
                        s = dist;
                        surfID = surf;
#ifdef MONITOR
                        std::cout << "\tThis distance is accepted as the new best estimate." << std::endl;
#endif
                    }
#ifdef MONITOR
                    else std::cout << "\tThis distance is rejected." << std::endl;
#endif
                }
            } else{
                // Check surfaces phi = phi_{k-1/2} and phi_{k+1/2}
                double phi = surf == 4 ? mesh->axis(2)[k] : mesh->axis(2)[k+1];

                double B = omega[1]*std::sin(phi)/b - omega[2]*std::cos(phi)/c;
                if (B == 0) continue;
                double C = ray.position()[1]*std::sin(phi)/b - ray.position()[2]*std::cos(phi)/c;
                double dist = -C/B;

                // Bound checking the root
                Eigen::Vector3d dest = mesh->coordinates(ray.position() + omega*dist);
#ifdef MONITOR
                std::cout << "On surface phi = " << phi << ", one possible distance is s = " << dist << ". ";
                std::cout << "The destination coordinates are " << dest.transpose() << std::endl;
#endif
                bool validS = dist > eps && dist < s;
                bool validR = dest[0] >= mesh->axis(0)[i] && dest[0] <= mesh->axis(0)[i+1];
                bool validMu = dest[1] >= mesh->axis(1)[j] && dest[1] <= mesh->axis(1)[j+1];
                bool validPhi = std::abs(dest[2]-phi) < mconst::pi; // make sure that it maps to the surface phi and not phi +/- pi
                if (validS && validR && validMu && validPhi){
                    s = dist;
                    surfID = surf;
#ifdef MONITOR
                    std::cout << "\tThis distance is accepted as the new best estimate." << std::endl;
#endif
                }
#ifdef MONITOR
                else std::cout << "\tThis distance is rejected." << std::endl;
#endif
            }
        }
        if (surfID == 6) throw std::runtime_error("ERROR: Ray tracing failure. Cannot solve for the next cell in the trajectory.");
#ifdef MONITOR
        else std::cout << "Ray will travel a distance of " << s << " and land on surface #" << surfID << std::endl;
#endif

        // Update heat source
        double alpha = _mat->computeProperty(Prop::attenuationCoefficient, vars);
        double Ed = ray.intensity()*(1 - std::exp(-alpha*s)); // deposited energy;
        double rho = _mat->computeProperty(Prop::density, vars);
        double Cp = _mat->computeProperty(Prop::heatCapacity, vars);
        double V = mesh->volume(i,j,k);
        q(i,j,k) = Ed/(rho*Cp*V);
        
        ray.setIntensity(ray.intensity() - Ed); // remaining energy
        if (ray.intensity() < IThreshold*_I0){
            if (dist(rng) < ray.intensity()/_I0) ray.setIntensity(IThreshold*_I0); // ray survives
            else{
                ray.setIntensity(0); // ray does not survive, heat source calculation is complete
                return q.matrix();
            }
        }

        // Advance the ray forward
        Eigen::Vector3d nextPosition = ray.position() += ray.direction().value()*s;
        ray.setPoistion(nextPosition);
        coord = mesh->coordinates(ray.position());
        n1 = n2;

        // Terminating condition
        normal = mesh->normalVector(0, coord[0], coord[1], coord[2]);
        if (surfID % 2 == 0) normal.value() *= -1; // lower surfaces, invert the normal
        if (surfID == 1 && i == Nr-1){
            // Ray arrived at an exit, deciding if it actually exits or is reflected back in
            n2 = 1.0;
            UnitVector3d transmit = ray.direction().refract(normal, n1, n2);
            if (transmit){
                ray.direction() = transmit;
#ifdef MONITOR2
                std::cout << "[" << ray.position().transpose() << "] on cell " << i << ", " << j << ", " << k << std::endl;
#endif
                return q.matrix();
            } else{
                // Total internal reflection, since there is no refraction
#ifdef MONITOR
                std::cout << "Total internal reflection: ray is reflected back to the medium." << std::endl;
                std::cout << "The normal vector is " << normal << std::endl;
#endif
                ray.setDirection(ray.direction().reflect(normal));
                omega = ray.direction().value();
            }
        } else{
            if (surfID == 0) i--;
            else if (surfID == 1) i++;
            else if (surfID == 2) j--;
            else if (surfID == 3) j++;
            else if (surfID == 4) k = k == 0 ? Nphi-1 : k-1;
            else k = k == Nphi-1 ? 0 : k+1;

            vars.at(_var) = u(i,j,k);
            n2 = _mat->computeProperty(Prop::refractiveIndex, vars);
            UnitVector3d transmit = ray.direction().refract(normal, n1, n2);
            if (transmit) ray.setDirection(ray.direction().refract(normal, n1, n2));
            /* Forced transmission if Snell's Law can't be solved for - i.e. zeroth-order approximation of the Eikonal equation */
            omega = ray.direction().value();
        }
#ifdef MONITOR
    std::cout << std::endl;
#endif
    }
}