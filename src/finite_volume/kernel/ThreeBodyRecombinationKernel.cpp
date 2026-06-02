#include "ThreeBodyRecombinationKernel.h"
#include "Constants.h"
#include "ElectronTemperatureAux.h"
#include "GaussLaguerre.h"
#include "Material.h"
#include "SpatialMesh.h"
#include "StateMesh.h"

Eigen::MatrixXd ThreeBodyRecombinationKernel::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the x axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the y axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the z axis
    Eigen::MatrixXd q(6, Nx*Ny*Nz);

    const double v = std::sqrt(2/pconst::m_e);
    double sahaFactor = 1.0; // // assuming g_bound / (2*g_ion) = 1.0 for now
    for (Eigen::Index i = 0; i < Nx; i++){
        for (Eigen::Index j = 0; j < Ny; j++){
            for (Eigen::Index k = 0; k < Nz; k++){
                auto vars = u.stateMap(i,j,k);
                double Eb = _mat->computeProperty("binding_energy", vars);
                double r = pconst::Ry / Eb;
                double uu = _mat->computeProperty("orbital_kinetic_energy", vars) / Eb;
                double kT = pconst::k_B*ElectronTemperatureAux().computeValue(vars);
                if (kT <= 0.0) throw std::invalid_argument("ERROR: Non-positive electron temperature encountered.");
                double pdf = 2/std::sqrt(mconst::pi) * std::pow(kT, -1.5);

                auto xsBEB = [&](const double& E){
                    double t = E/Eb;
                    double A = 4*mconst::pi * (pconst::a0*pconst::a0) * r*r;
                    double B = t + uu + 1;
                    double C = std::log(t)*(1.0-1/t/t)/2 + 1 - 1/t - std::log(t)/(t+1);
                    return A/B*C;
                };

                auto func = [&](const double& E){ return xsBEB(E) * v*pdf*E; };

                GaussLaguerre quad(8, 1.0/kT); // Use ungeneralized Gauss-Laguerre because the cross section function may have a discontinuity at E = 0
                Quadrature<>::IntegrationDomain D{Eb, mconst::infty}; // lower bound is Eb since that's the threshold of ionization
                double S_ci = quad.integrate(func, D);

                double lambda_th = pconst::h / std::sqrt(2*mconst::pi * pconst::m_e * kT);
                double K = 1.0/(lambda_th*lambda_th*lambda_th) * sahaFactor * std::exp(-Eb/kT);
                double S_3br = S_ci/K;

                double ne = vars.at("electron_number_density");
                double ni = vars.at("ion_number_density");
                double Ei = vars.at("ion_energy")/ni;
                double Gamma_3br = S_3br*ne*ne*ni;

                Eigen::Index ind = (i*Nx + j)*Ny + k;
                q(0, ind) = 0.0; // Gamma_3br * (_mat->computeProperty("molecular_mass") / pconst::N_A); // neutral mass density
                q(1, ind) = 0.0; // Gamma_3br*Ei; // neutral energy
                q(2, ind) = -Gamma_3br; // electron number density
                q(3, ind) = Gamma_3br*Eb; // electron energy
                q(4, ind) = -Gamma_3br; // ion number density
                q(5, ind) = -Gamma_3br*Ei; // ion energy
            }
        }
    }
    return q;
}

Eigen::VectorXi ThreeBodyRecombinationKernel::stateID(const StateMesh& u) const noexcept{
    Eigen::VectorXi ind(6);
    ind(0) = u.stateID("density");
    ind(1) = u.stateID("energy");
    ind(2) = u.stateID("electron_number_density");
    ind(3) = u.stateID("electron_energy");
    ind(4) = u.stateID("ion_number_density");
    ind(5) = u.stateID("ion_energy");
    return ind;
}