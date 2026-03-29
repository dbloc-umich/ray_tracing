#include "ElectronHPEnergyTransferKernel.h"
#include "Constants.h"
#include "EquationOfState.h"
#include "SpatialMesh.h"
#include "StateMesh.h"
#include "TemperatureFromEnergyAux.h"

ElectronHPEnergyTransferKernel::ElectronHPEnergyTransferKernel(std::shared_ptr<Material> mat, std::shared_ptr<EquationOfState> eos):
    Kernel(mat),
    _eos(eos)
{
    assert(_eos);
}

Eigen::MatrixXd ElectronHPEnergyTransferKernel::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the r axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the mu axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the phi axis
    //Eigen::Index Ns   = _s.size();           // number of variables that this residual affects, should be 2 for now
    Eigen::MatrixXd q(2, Nx*Ny*Nz);

    TemperatureFromEnergyAux ThAux(_eos);
    for (Eigen::Index i = 0; i < Nx; i++){
        for (Eigen::Index j = 0; j < Ny; j++){
            for (Eigen::Index k = 0; k < Nz; k++){
                auto vars = u.matProp(i,j,k);
                double Ee = vars.at("electron_energy");
                double f = (2*pconst::m_e*pconst::N_A)/_mat->computeProperty("molecular_mass", vars); // conversion factor from momentum to energy relaxation frequency
                double nu_ei = _mat->computeProperty("electron_ion_collision_frequency", vars)*f;
                double nu_en = _mat->computeProperty("electron_neutral_collision_frequency", vars)*f;
                double ne = _mat->computeProperty("electron_density", vars);
                double Th = ThAux.computeValue(vars);
                double C = 1.5*ne*pconst::k_B;
                double qeh = (nu_ei+nu_en)*(Ee-C*Th);
                // std::cout << "Ee = " << Ee << ", T = " << Th << ", E = " << C*Th << std::endl;
                q(0, (i*Ny + j)*Nz + k) = qeh; // neutral energy
                q(1, (i*Ny + j)*Nz + k) = -qeh; // electron energy
            }
        }
    }
    return q;
}

Eigen::VectorXi ElectronHPEnergyTransferKernel::stateID(const StateMesh& u) const noexcept{
    Eigen::VectorXi ind(2);
    ind(0) = u.stateID("energy");
    ind(1) = u.stateID("electron_energy");
    return ind;
}