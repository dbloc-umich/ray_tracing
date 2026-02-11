#include "FVLaserSource.h"
#include "Constants.h"
#include "Ellipsoid.h"
#include "FVSpatialMesh.h"
#include "FVStateMesh.h"
#include "Ray.h"
#include "Sphere.h"

static constexpr double eps = 1.0e-9;

FVLaserSource::FVLaserSource(std::shared_ptr<Material> mat, Ray& ray, Shape* shape):
    FVResidual(mat, PropVariable::temperature),
    _ray(ray),
    _shape(shape)
{
    assert(!_shape->surfaceContains(_ray.position()) && !_shape->encloses(_ray.position()));
}

Eigen::ArrayXd FVLaserSource::computeResidual(const FVStateMesh& u) const{
    auto mesh = u.meshPtr();
    Eigen::Index Nr   = mesh->axisSize(0)-1; // number of cells on the r axis
    Eigen::Index Nmu  = mesh->axisSize(1)-1; // number of cells on the mu axis
    Eigen::Index Nphi = mesh->axisSize(2)-1; // number of cells on the phi axis
    FVStateMesh q(mesh, 0);

    double a, b, c; // semi-axes of the ellipse
    if (Sphere* sph = dynamic_cast<Sphere*>(_shape)){
        a = sph->radius();
        b = sph->radius();
        c = sph->radius();
    } else if (Ellipsoid* ell = dynamic_cast<Ellipsoid*>(_shape)){
        a = std::max({ell->semiAxis(0), ell->semiAxis(1), ell->semiAxis(2)});
        c = std::min({ell->semiAxis(0), ell->semiAxis(1), ell->semiAxis(2)});
        b = ell->semiAxis(0) + ell->semiAxis(1) + ell->semiAxis(2) - a - c;
    } else throw std::invalid_argument("ERROR: Laser heat source has only been implemented for spherical and ellipsoidal meshes");

    Ray ray = _ray;
    double s = _shape->distanceToSurface(ray.position(), ray.direction());
    if (std::isnan(s)) return q.array(); // ray does not enter shape

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
    UnitVector3d normal = mesh->normalVector(0, ray.position()[0], ray.position()[1], ray.position()[2]);
    double n1 = 1.0;
    Material::PropVars vars{{_var, u(i,j,k)}};
    double n2 = _mat->computeProperty(Prop::refractiveIndex, vars);
    ray.setDirection(ray.direction().refract(normal, n1, n2));
    omega = ray.direction().value();

    bool inside = true;
    while (inside){
        s = std::numeric_limits<double>::max();
        std::size_t surfID = 6; // surfID of the next surface, 6 is an invalid ID

        // Find the exit surface and its distance
        for (std::size_t surf = 0; surf < 6; surf++){
            if (surf == 0 || surf == 1){
                // Check surfaces r = r_{i-1/2} and r_{i+1/2}
                double r = surf == 0 ? mesh->axis(0)[i] : mesh->axis(0)[i+1];
                if (surf == 0 && i == 0) continue; // can theoretically passes through the origin, but it's computationally impossible
                Eigen::DiagonalMatrix<double, 3> M{1/(a*a), 1/(b*b), 1/(c*c)};
                double A = omega.dot(M*omega);
                double B = omega.dot(M*ray.position());
                double C = ray.position().dot(M*ray.position()) - r*r;
                
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
                    bool validS = dist > eps && dist < s;
                    bool validMu = dest[1] >= mesh->axis(1)[j] && dest[1] <= mesh->axis(1)[j+1];
                    bool validPhi = dest[2] >= mesh->axis(2)[k] && dest[2] <= mesh->axis(2)[k+1];
                    if (validS && validMu && validPhi){
                        s = dist;
                        surfID = surf;
                    }
                }
            } else if (surf == 2 || surf == 3){
                // Check surfaces mu = mu_{j-1/2} and mu_{j+1/2}
                double mu = surf == 2 ? mesh->axis(1)[j] : mesh->axis(1)[j+1];
                if ((surf == 2 && j == 0) || (surf == 3 && j == Nmu-1)) continue; // can theoretically passes through the major axis, but it's computationally impossible

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
                    bool validS = dist > eps && dist < s;
                    bool validR = dest[0] >= mesh->axis(0)[i] && dest[0] <= mesh->axis(0)[i+1];
                    bool validMu = dest[0]*mu >= 0; // make sure that it maps to the surface mu and not -mu
                    bool validPhi = dest[2] >= mesh->axis(2)[k] && dest[2] <= mesh->axis(2)[k+1];
                    if (validS && validR && validMu && validPhi){
                        s = dist;
                        surfID = surf;
                    }
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
                bool validS = dist > eps && dist < s;
                bool validR = dest[0] >= mesh->axis(0)[i] && dest[0] <= mesh->axis(0)[i+1];
                bool validMu = dest[1] >= mesh->axis(1)[j] && dest[1] <= mesh->axis(1)[j+1];
                bool validPhi = std::abs(dest[2]-phi) < mconst::pi; // make sure that it maps to the surface phi and not phi +/- pi
                if (validS && validR && validMu && validPhi){
                    s = dist;
                    surfID = surf;
                }
            }
        }
        if (surfID == 6) throw std::runtime_error("ERROR: Ray tracing failure. Cannot solve for the next cell in the trajectory.");
        
        // Update heat source
        double alpha = _mat->computeProperty(Prop::attenuationCoefficient, vars);
        double Ed = ray.intensity()*(1 - std::exp(-alpha*s)); // deposited energy;
        double rho = _mat->computeProperty(Prop::density, vars);
        double Cp = _mat->computeProperty(Prop::heatCapacity, vars);
        double V = mesh->volume(i,j,k);
        q(i,j,k) = Ed/(rho*Cp*V);
        ray.setIntensity(ray.intensity() - Ed); // remaining energy
        
        // Advance the ray forward
        Eigen::Vector3d nextPosition = ray.position() += ray.direction().value()*s;
        ray.setPoistion(nextPosition);
        n1 = n2;

        // Terminating condition
        if (surfID == 1 && i == Nr-1){
            n2 = 1.0;
            normal =  mesh->normalVector(0, ray.position()[0], ray.position()[1], ray.position()[2]);
            ray.setDirection(ray.direction().refract(normal, n1, n2));
            inside = false;
        } else{
            Eigen::Vector3d coord = mesh->coordinates(ray.position());
            if (surfID == 0){
                i--;
                normal = -mesh->normalVector(0, coord[0], coord[1], coord[2]);
            } else if (surfID == 1){
                i++;
                normal = mesh->normalVector(0, coord[0], coord[1], coord[2]);
            } else if (surfID == 2){
                j--;
                normal = -mesh->normalVector(1, coord[0], coord[1], coord[2]);
            } else if (surfID == 3){
                j++;
                normal = mesh->normalVector(1, coord[0], coord[1], coord[2]);
            } else if (surfID == 4){
                k = k == 0 ? Nphi-1 : k-1;
                normal = -mesh->normalVector(2, coord[0], coord[1], coord[2]);
            } else{
                k = k == Nphi-1 ? 0 : k+1;
                normal = mesh->normalVector(2, coord[0], coord[1], coord[2]);
            }

            vars.at(_var) = u(i,j,k);
            double n2 = _mat->computeProperty(Prop::refractiveIndex, vars);
            ray.setDirection(ray.direction().refract(normal, n1, n2));
            omega = ray.direction().value();
        }
    }
    return q.array();
}