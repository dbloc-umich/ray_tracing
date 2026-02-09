#include "FVLaserSource.h"
#include "FVSpatialMesh.h"
#include "FVStateMesh.h"
#include "NewtonSolver.h"
#include "Ray.h"
#include "Shape.h"

FVLaserSource::FVLaserSource(std::shared_ptr<Material> mat, Ray& ray, Shape* shape):
    FVResidual(mat, PropVariable::temperature),
    _ray(ray),
    _shape(shape)
{
    assert(!_shape->surfaceContains(_ray.position()) && !_shape->encloses(_ray.position()));
}

FVStateMesh FVLaserSource::computeResidual(const FVStateMesh& u) const{
    auto mesh = u.meshPtr();
    FVStateMesh q(mesh);

    double s = _shape->distanceToSurface(_ray.position(), _ray.direction());
    if (std::isnan(s)) return q; // ray does not enter shape

    Eigen::Vector3d omega = _ray.direction().value();
    _ray.setPoistion(_ray.position() + omega*s);

    /*
    convert current ray position into the shape coordinates   
    */
    
    // Locates the first cell that the ray enters
    Eigen::Index Nr = mesh->axis(0).size()-1; // number of cells on the r axis
    Eigen::Index Nmu = mesh->axis(1).size()-1; // number of cells on the mu axis
    Eigen::Index Nphi = mesh->axis(2).size()-1; // number of cells on the phi axis

    Eigen::Vector3d coord = mesh->coordinates(_ray.position()[0], _ray.position()[1], _ray.position()[2]);
    Eigen::Index i = Nr-1; // only works for spherical and ellipsoidal coordinates
    Eigen::Index j = std::lower_bound(mesh->axis(1).cbegin(), mesh->axis(1).cend(), coord[1]) - mesh->axis(1).cbegin() - 1;
    Eigen::Index k = std::lower_bound(mesh->axis(2).cbegin(), mesh->axis(2).cend(), coord[2]) - mesh->axis(2).cbegin() - 1;

    // Refract the ray
    UnitVector3d normal = mesh->normalVector(0, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
    double n1 = 1.0;
    Material::PropVars vars{{PropVariable::temperature, u(i,j,k)}};
    double n2 = _mat->computeProperty(Prop::refractiveIndex, vars);
    _ray.setDirection(_ray.direction().refract(normal, n1, n2));
    omega = _ray.direction().value();

    bool inside = true;
    NewtonSolver<3> newton;
    std::size_t onSurf = 1;
    while (inside){
        s = std::numeric_limits<double>::max();
        std::size_t surfID = 6; // surfID of the next surface, 6 is an invalid ID

        // Find the exit surface and its distance
        for (std::size_t surf = 0; surf < 6; surf++){
            if (surf == onSurf) continue;
            if (surf == 0 || surf == 1){
                // Check surfaces r = r_{i-1/2} and r_{i+1/2}
                double r = surf == 0 ? mesh->axis(0)[i] : mesh->axis(0)[i+1];
                if (surf == 0 && i == 0) continue; // can theoretically passes through the origin, but it's computationally impossible
                
                double muL = mesh->axis(1)[j];
                double muU = mesh->axis(1)[j+1];
                double phiL = mesh->axis(2)[k];
                double phiU = mesh->axis(2)[k+1];
                auto f = [&, this](const Eigen::Vector3d& v) -> Eigen::Vector3d
                {
                    double dist = v[0];
                    double muOut = v[1];
                    double phiOut = v[2];
                    return _ray.position() + omega*dist - mesh->cartesian(r, muOut, phiOut);
                };
                auto J = [&, this](const Eigen::Vector3d& v){
                    Eigen::Matrix3d jac;
                    double muOut = v[1];
                    double phiOut = v[2];
                    jac.col(0) = omega;
                    jac.col(1) = -mesh->basisVector(1, r, muOut, phiOut);
                    jac.col(2) = -mesh->basisVector(2, r, muOut, phiOut);
                    return jac;
                };
                Eigen::Vector3d x{mesh->axis(0)[i+1]-mesh->axis(0)[i], (muL+muU)/2, (phiL+phiU)/2};
                newton.setFunction(f, J);
                auto status = newton.solve(x);
                if (status == NLStatus::Success){
                    bool validS = x[0] > 0 && x[0] < s;
                    bool validMu = x[1] >= muL && x[1] <= muU;
                    bool validPhi = x[2] >= phiL && x[2] <= phiU;
                    if (validS && validMu && validPhi){
                        s = x[0];
                        surfID = surf;
                    }
                }
            } else if (surf == 2 || surf == 3){
                // Check surfaces mu = mu_{j-1/2} and mu_{j+1/2}
                double mu = surf == 2 ? mesh->axis(1)[j] : mesh->axis(1)[j+1];
                if ((surf == 2 && j == 0) || (surf == 3 && j == Nmu-1)) continue; // can theoretically passes through the major axis, but it's computationally impossible

                double rL = mesh->axis(0)[i];
                double rU = mesh->axis(0)[i+1];
                double phiL = mesh->axis(2)[k];
                double phiU = mesh->axis(2)[k+1];
                auto f = [&, this](const Eigen::Vector3d& v) -> Eigen::Vector3d {
                    double dist = v[0];
                    double rOut = v[1];
                    double phiOut = v[2];
                    return _ray.position() + omega*dist - mesh->cartesian(rOut, mu, phiOut);
                };
                auto J = [&, this](const Eigen::Vector3d& v){
                    Eigen::Matrix3d jac;
                    double rOut = v[1];
                    double phiOut = v[2];
                    jac.col(0) = omega;
                    jac.col(1) = -mesh->basisVector(0, rOut, mu, phiOut);
                    jac.col(2) = -mesh->basisVector(2, rOut, mu, phiOut);
                    return jac;
                };
                Eigen::Vector3d x{mesh->axis(1)[j+1]-mesh->axis(1)[j], (rL+rU)/2, (phiL+phiU)/2};
                newton.setFunction(f, J);
                auto status = newton.solve(x);
                if (status == NLStatus::Success){
                    bool validS = x[0] > 0 && x[0] < s;
                    bool validR = x[1] >= rL && x[1] <= rU;
                    bool validPhi = x[2] >= phiL && x[2] <= phiU;
                    if (validS && validR && validPhi){
                        s = x[0];
                        surfID = surf;
                    }
                }
            }
            else{
                // Check surfaces phi = phi_{k-1/2} and phi_{k+1/2}
                double phi = surf == 4 ? mesh->axis(2)[k] : mesh->axis(2)[k+1];

                double rL = mesh->axis(0)[i];
                double rU = mesh->axis(0)[i+1];
                double muL = mesh->axis(1)[j];
                double muU = mesh->axis(1)[j+1];
                auto f = [&, this](const Eigen::Vector3d& v) -> Eigen::Vector3d
                {
                    double dist = v[0];
                    double rOut = v[1];
                    double muOut = v[2];
                    return _ray.position() + omega*dist - mesh->cartesian(rOut, muOut, phi);
                };
                auto J = [&, this](const Eigen::Vector3d& v){
                    Eigen::Matrix3d jac;
                    double rOut = v[1];
                    double muOut = v[2];
                    jac.col(0) = omega;
                    jac.col(1) = -mesh->basisVector(0, rOut, muOut, phi);
                    jac.col(2) = -mesh->basisVector(1, rOut, muOut, phi);
                    return jac;
                };
                Eigen::Vector3d x{mesh->axis(2)[k+1]-mesh->axis(2)[k], (rL+rU)/2, (muL+muU)/2};
                newton.setFunction(f, J);
                auto status = newton.solve(x);
                if (status == NLStatus::Success){
                    bool validS = x[0] > 0 && x[0] < s;
                    bool validR = x[1] >= rL && x[1] <= rU;
                    bool validMu = x[2] >= muL && x[2] <= muU;
                    if (validS && validR && validMu){
                        s = x[0];
                        surfID = surf;
                    }
                }
            }
        }
        if (surfID == 6) throw std::runtime_error("ERROR: Ray tracing failure. Cannot solve for the next cell in the trajectory.");

        // Update heat source
        double alpha = _mat->computeProperty(Prop::attenuationCoefficient, vars);
        double Ed = _ray.intensity()*(1 - std::exp(-alpha*s)); // deposited energy;
        double rho = _mat->rho(101325, u(i,j,k));
        double Cp = _mat->Cp(101325, u(i,j,k));
        q(i,j,k) = Ed/(rho*Cp);
        _ray.setIntensity(_ray.intensity() - Ed); // remaining energy
        
        // Advance the ray forward
        Eigen::Vector3d nextPosition = _ray.position() += _ray.direction().value()*s;
        _ray.setPoistion(nextPosition);
        n1 = n2;

        // Terminating condition
        if (surfID == 1 && i == Nr-1){
            n2 = 1.0;
            normal =  mesh->normalVector(0, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
            _ray.setDirection(_ray.direction().refract(normal, n1, n2));
            inside = false;
        } else{
            if (surfID == 0){
                i--;
                normal = -mesh->normalVector(0, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
                onSurf = 1;
            } else if (surfID == 1){
                i++;
                normal = mesh->normalVector(0, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
                onSurf = 0;
            } else if (surfID == 2){
                j--;
                normal = -mesh->normalVector(1, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
                onSurf = 3;
            } else if (surfID == 3){
                j++;
                normal = mesh->normalVector(1, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
                onSurf = 2;
            } else if (surfID == 4){
                k = k == 0 ? Nphi-1 : k-1;
                normal = -mesh->normalVector(2, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
                onSurf = 5;
            } else{
                k = k == Nphi-1 ? 0 : k+1;
                normal = mesh->normalVector(2, _ray.position()[0], _ray.position()[1], _ray.position()[2]);
                onSurf = 4;
            }

            vars.at(PropVariable::temperature) = u(i,j,k);
            double n2 = _mat->computeProperty(Prop::refractiveIndex, vars);
            _ray.setDirection(_ray.direction().refract(normal, n1, n2));
            omega = _ray.direction().value();
        }
    }
    return q;
}