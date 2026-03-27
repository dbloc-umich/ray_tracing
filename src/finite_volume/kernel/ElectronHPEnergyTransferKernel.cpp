#include "ElectronHPEnergyTransferKernel.h"
#include "Constants.h"
#include "ElectronTemperatureAux.h"
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
    Eigen::MatrixXd q = Eigen::MatrixXd::Zero(2, Nx*Ny*Nz);

    TemperatureFromEnergyAux ThAux(_eos);
    ElectronTemperatureAux TeAux;

    for (Eigen::Index i = 0; i < Nx; i++){
        for (Eigen::Index j = 0; j < Ny; j++){
            for (Eigen::Index k = 0; k < Nz; k++){
                auto vars = u.matProp(i,j,k);
                double Te = TeAux.computeValue(vars);
                if (!std::isnan(Te)){
                    double nu_ei = _mat->computeProperty("electron_ion_collision_frequency", vars);
                    double nu_en = _mat->computeProperty("electron_neutral_collision_frequency", vars);
                    double ne = _mat->computeProperty("electron_density", vars);
                    double Th = ThAux.computeValue(vars);
                    double qeh = 1.5*(nu_ei+nu_en)*ne*pconst::k_B*(Te-Th);
                    // std::cout << "Te-Th = " << Te-Th << ", nu = " << nu_ei+nu_en << ", ne = " << ne << ", qeh = " << qeh << std::endl;
                    q(0, (i*Ny + j)*Nz + k) = qeh;
                    q(1, (i*Ny + j)*Nz + k) = -q(0, (i*Ny + j)*Nz + k);
                }
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