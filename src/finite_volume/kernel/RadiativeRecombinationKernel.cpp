#include "RadiativeRecombinationKernel.h"
#include "Constants.h"
#include "ElectronTemperatureAux.h"
#include "GaussLaguerre.h"
#include "GeneralizedGaussLaguerre1.h"
#include "Material.h"
#include "SpatialMesh.h"
#include "StateMesh.h"

Eigen::MatrixXd RadiativeRecombinationKernel::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the x axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the y axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the z axis
    Eigen::MatrixXd q(6, Nx*Ny*Nz);

    double milneFactor = 1.0; // assuming g_bound / (2*g_ion) = 1.0 for now
    double v = std::sqrt(2/pconst::m_e);
    for (Eigen::Index i = 0; i < Nx; i++){
        for (Eigen::Index j = 0; j < Ny; j++){
            for (Eigen::Index k = 0; k < Nz; k++){
                auto vars = u.stateMap(i,j,k);
                double Eb = _mat->computeProperty("binding_energy", vars);
                double kT = pconst::k_B*ElectronTemperatureAux().computeValue(vars);
                if (kT <= 0.0) throw std::invalid_argument("ERROR: Non-positive electron temperature encountered.");
                double pdf = 2/std::sqrt(mconst::pi) * std::pow(kT, -1.5);

                auto xsPi = [&](const double& E){
                    double nu = (E + Eb)/pconst::h; // photon frequency in s-1
                    double lambda = pconst::c / nu; // photon vacuum wavelength in m
                    vars["wavelength"] = lambda;
                    return _mat->computeProperty("photoionization_cross_section", vars);
                };

                // The integral of S_rr but without the factor of exp(-E/kT), which is the implicit weight of ungeneralized Gauss-Laguerre
                auto func = [&](const double& E){
                    double hnu = E + Eb;
                    double mc2 = pconst::m_e * pconst::c * pconst::c;
                    double C = milneFactor * hnu*hnu/mc2;
                    return xsPi(E) * v*pdf*C;
                };

                Quadrature<>::IntegrationDomain D{0.0, mconst::infty};
                // Reaction integral
                GaussLaguerre GL0(8, 1.0/kT);
                double S_rr = GL0.integrate(func, D);

                // Energy integral
                GeneralizedGaussLaguerre1 GL1(8, 1.0/kT);
                double K_rr = GL1.integrate(func, D);

                double ne = vars.at("electron_number_density");
                double ni = vars.at("ion_number_density");
                double Ei = vars.at("ion_energy")/ni;
                double Gamma_rr = S_rr*ne*ni;
                double Qe_rr = K_rr*ne*ni;
                // std::cout << "Gamma_rr = " << Gamma_rr << ", Qe_rr = " << Qe_rr << std::endl;

                Eigen::Index ind = (i*Nx + j)*Ny + k;
                q(0, ind) = 0.0; // Gamma_ci * (_mat->computeProperty("molecular_mass") / pconst::N_A); // neutral mass density
                q(1, ind) = 0.0; // Gamma_ci*Ei; // neutral energy
                q(2, ind) = -Gamma_rr; // electron number density
                q(3, ind) = -Qe_rr; // electron energy
                q(4, ind) = -Gamma_rr; // ion number density
                q(5, ind) = -Gamma_rr*Ei; // ion energy
            }
        }
    }
    return q;
}

Eigen::VectorXi RadiativeRecombinationKernel::stateID(const StateMesh& u) const noexcept{
    Eigen::VectorXi ind(6);
    ind(0) = u.stateID("density");
    ind(1) = u.stateID("energy");
    ind(2) = u.stateID("electron_number_density");
    ind(3) = u.stateID("electron_energy");
    ind(4) = u.stateID("ion_number_density");
    ind(5) = u.stateID("ion_energy");
    return ind;
}