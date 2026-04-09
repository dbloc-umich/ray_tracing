#include "CollisionalIonizationKernel.h"
#include "Constants.h"
#include "ElectronTemperatureAux.h"
#include "GaussLaguerre.h"
#include "SpatialMesh.h"
#include "StateMesh.h"

Eigen::MatrixXd CollisionalIonizationKernel::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the r axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the mu axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the phi axis
    Eigen::MatrixXd q(3, Nx*Ny*Nz);

    for (Eigen::Index i = 0; i < Nx; i++){
        for (Eigen::Index j = 0; j < Ny; j++){
            for (Eigen::Index k = 0; k < Nz; k++){
                auto vars = u.stateMap(i,j,k);
                double Eb = _mat->computeProperty("binding_energy", vars);
                double r = pconst::Ry / Eb;
                double uu = _mat->computeProperty("orbital_kinetic_energy", vars) / Eb;
                double kT = pconst::k_B*ElectronTemperatureAux().computeValue(vars);

                auto xsBEB = [&](const double& E){
                    double t = E/Eb;
                    double A = 4*mconst::pi * (pconst::a0*pconst::a0) * r*r;
                    double B = t + uu + 1;
                    double C = std::log(t)*(1.0-1/t/t)/2 + 1 - 1/t - std::log(t)/(t+1);
                    return A/B*C;
                };

                auto func = [&](const double& E){
                    double v = std::sqrt(2/pconst::m_e);
                    double pdf = 2/std::sqrt(mconst::pi) * std::pow(kT, -1.5);
                    return xsBEB(E) * v*pdf*E;
                };

                GaussLaguerre quad(8, 1.0/kT); // Use ungeneralized Gauss-Laguerre because the cross section function may be a discontinuity at E = 0
                Quadrature<>::IntegrationDomain D{Eb, mconst::infty};
                double S_ci = quad.integrate(func, D);
                // std::cout << "S_ci = " << S_ci << ", Eb = " << Eb << ", u = " << u << ", Te = " << kT/pconst::k_B << std::endl;
                double ne = _mat->computeProperty("electron_density", vars);
                double n = _mat->computeProperty("number_density", vars);
                double En = vars.at("energy")/n; // energy per neutral particle
                double Gamma_ci = S_ci*ne*n;

                Eigen::Index ind = (i*Nx + j)*Ny + k;
                q(0, ind) = -Gamma_ci*En;
                q(1, ind) = -Gamma_ci*Eb;
                q(2, ind) = Gamma_ci;
            }
        }
    }
    return q;
}

Eigen::VectorXi CollisionalIonizationKernel::stateID(const StateMesh& u) const noexcept{
    Eigen::VectorXi ind(3);
    ind(0) = u.stateID("energy");
    ind(1) = u.stateID("electron_energy");
    ind(2) = u.stateID("ion_number_density");
    return ind;
}