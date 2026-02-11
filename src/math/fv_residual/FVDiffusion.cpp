#include "FVDiffusion.h"
#include "FVBoundaryCondition.h"
#include "FVSpatialMesh.h"
#include "FVStateMesh.h"

FVDiffusion::FVDiffusion(const std::vector<std::shared_ptr<FVBoundaryCondition>>& bc, std::shared_ptr<Material> mat, PropVariable var, Prop prop):
    FVResidual(mat, var),
    _bc(bc),
    _prop(prop)
{
    assert(_bc.size() == 6);
    assert(var == PropVariable::temperature || var == PropVariable::velocity);
    for (std::size_t i = 0; i < 3; i++){
        FVBoundaryCondition* pL = dynamic_cast<FVBoundaryCondition*>(_bc[2*i].get());
        FVBoundaryCondition* pU = dynamic_cast<FVBoundaryCondition*>(_bc[2*i+1].get());
        if (pL && pU) _periodic.push_back(true);
        else if (!pL && !pU) _periodic.push_back(false);
        else throw std::invalid_argument("ERROR: The two boundaries of a variable have to be either both periodic or both non-periodic.");
    }
}

Eigen::ArrayXd FVDiffusion::computeResidual(const FVStateMesh& u) const{
    auto mesh = u.meshPtr();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the 1st axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the 2nd axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the 3rd axis
    FVStateMesh F(mesh);

    for (Eigen::Index i = 0; i < Nx; i++){
        for (Eigen::Index j = 0; j < Ny; j++){
            for (Eigen::Index k = 0; k < Nz; k++){
                double u0 = u(i,j,k);
                double rho = _mat->computeProperty(Prop::density, {{_var, u0}});
                double Cp = _mat->computeProperty(Prop::heatCapacity, {{_var, u0}});
                double V = mesh->volume(i,j,k);
                double D = (_prop == Prop::none) ? 1.0 : _mat->computeProperty(_prop, {{_var, u0}});

                // Flux from and to neighboring cells
                std::vector<Eigen::Index> surfaces;
                if (i < Nx-1 || (i == Nx-1 && _periodic[0])) surfaces.push_back(1);
                if (j < Ny-1 || (j == Ny-1 && _periodic[1])) surfaces.push_back(3);
                if (k < Nz-1 || (k == Nz-1 && _periodic[2])) surfaces.push_back(5);
                for (auto surfID: surfaces){
                    // Only considering the "upper" surfaces, because of symmetry that the lower surfaces can be updated on a different cell
                    double uk, dr, Vk;
                    if (surfID == 1){
                        if (i < Nx-1){ // interior
                            uk = u(i+1,j,k);
                            dr = (mesh->axis(0)[i+2] - mesh->axis(0)[i])/2;
                            Vk = mesh->volume(i+1,j,k);
                        } else{ // periodic
                            uk = u(0,j,k);
                            dr = (mesh->axis(0)[i+1] - mesh->axis(0)[i])/2 + (mesh->axis(0)[1] - mesh->axis(0)[0])/2;
                            Vk = mesh->volume(0,j,k);
                        }
                    } else if (surfID == 3){
                        if (j < Ny-1){ // interior
                            uk = u(i,j+1,k);
                            dr = (mesh->axis(1)[j+2] - mesh->axis(1)[j])/2;
                            Vk = mesh->volume(i,j+1,k);
                        } else{ // periodic
                            uk = u(i,0,k);
                            dr = (mesh->axis(1)[j+1] - mesh->axis(1)[j])/2 + (mesh->axis(1)[1] - mesh->axis(1)[0])/2;
                            Vk = mesh->volume(i,0,k);
                        }
                    } else{
                        if (k < Nz-1){ // interior
                            uk = u(i,j,k+1);
                            dr = (mesh->axis(2)[k+2] - mesh->axis(2)[k])/2;
                            Vk = mesh->volume(i,j,k+1);
                        } else{ // periodic
                            uk = u(i,j,0);
                            dr = (mesh->axis(2)[k+1] - mesh->axis(2)[k])/2 + (mesh->axis(2)[1] - mesh->axis(2)[0])/2;
                            Vk = mesh->volume(i,j,0);
                        }
                    }

                    double Dk = (_prop == Prop::none) ? 1.0 : _mat->computeProperty(_prop, {{_var, uk}});
                    double du = Dk*uk - D*u0;
                    double flux = mesh->gradientDotN(du/dr, i, j, k, surfID) * mesh->area(i, j, k, surfID);
                    F(i,j,k) += flux/(rho*Cp*V);

                    // Update the neighboring cell with flux through its "lower" surface
                    double rhok =  _mat->computeProperty(Prop::density, {{_var, uk}});
                    double Cpk = _mat->computeProperty(Prop::heatCapacity, {{_var, uk}});
                    if (surfID == 1) F(i < Nx-1 ? i+1 : 0, j, k) -= flux/(rhok*Cpk*Vk);
                    else if (surfID == 3) F(i, j < Ny-1 ? j+1 : 0, k) -= flux/(rhok*Cpk*Vk);
                    else F(i, j, k < Nz-1 ? k+1 : 0) -= flux/(rhok*Cpk*Vk);
                }

                // Flux from boundaries
                if (i == 0){
                    Eigen::Vector3d r{ mesh->axis(0)[i], (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                    F(i,j,k) += _bc[0]->computeFlux(u0, r) * mesh->area(i,j,k,0) / (rho*Cp*V);
                }
                if (i == Nx-1){
                    Eigen::Vector3d r{ mesh->axis(0)[i], (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                    F(i,j,k) += _bc[1]->computeFlux(u0, r) * mesh->area(i,j,k,1) / (rho*Cp*V);
                }
                if (j == 0){
                    Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, mesh->axis(1)[j], (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                    F(i,j,k) += _bc[2]->computeFlux(u0, r) * mesh->area(i,j,k,2) / (rho*Cp*V);
                }
                if (j == Ny-1){
                    Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, mesh->axis(1)[j], (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                    F(i,j,k) += _bc[3]->computeFlux(u0, r) * mesh->area(i,j,k,3) / (rho*Cp*V);
                }
                if (k == 0){
                    Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, mesh->axis(2)[k]};
                    F(i,j,k) += _bc[4]->computeFlux(u0, r) * mesh->area(i,j,k,4) / (rho*Cp*V);
                }
                if (k == Nz-1){
                    Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, mesh->axis(2)[k]};
                    F(i,j,k) += _bc[5]->computeFlux(u0, r) * mesh->area(i,j,k,5) / (rho*Cp*V);
                }
            }
        }
    }

    return F.array();
}