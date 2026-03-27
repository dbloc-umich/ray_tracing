#include "DiffusionKernel.h"
#include "AuxKernel.h"
#include "BoundaryCondition.h"
#include "SpatialMesh.h"
#include "StateMesh.h"
#include <iostream>

DiffusionKernel::DiffusionKernel(std::shared_ptr<Material> mat, std::string var,
                                 std::shared_ptr<AuxKernel> converter, std::string primitive, std::string prop):
    Kernel(mat),
    _var(std::move(var)),
    _converter(converter),
    _primitive(primitive),
    _prop(std::move(prop))
{}

Eigen::MatrixXd DiffusionKernel::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the 1st axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the 2nd axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the 3rd axis
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(1, Nx*Ny*Nz);
    Eigen::Index s = u.stateID(_var); // varID

    if (_prop != "none"){
        for (Eigen::Index i = 0; i < Nx; i++){
            for (Eigen::Index j = 0; j < Ny; j++){
                for (Eigen::Index k = 0; k < Nz; k++){
                    auto vars = u.matProp(i, j, k);
                    double u0 = u(s,i,j,k);
                    double T0 = _converter ? u0 : _converter->computeValue(vars);
                    vars[_primitive] = T0;
                    double D = _mat->computeProperty(_prop, vars);

                    // Flux from and to neighboring cells
                    std::vector<Eigen::Index> surfaces;
                    if (i < Nx-1 || (i == Nx-1 && u.bc(s,1)->isPeriodic())) surfaces.push_back(1);
                    if (j < Ny-1 || (j == Ny-1 && u.bc(s,3)->isPeriodic())) surfaces.push_back(3);
                    if (k < Nz-1 || (k == Nz-1 && u.bc(s,5)->isPeriodic())) surfaces.push_back(5);
                    for (auto surfID: surfaces){
                        // Only considering the "upper" surfaces, because of symmetry that the lower surfaces can be updated on a different cell
                        Eigen::Index in=i, jn=j, kn=k;
                        double dr;
                        if (surfID == 1){
                            if (i < Nx-1){ // interior
                                in = i+1;
                                dr = (mesh->axis(0)[i+2] - mesh->axis(0)[i])/2;
                            } else{ // periodic
                                in = 0;
                                dr = (mesh->axis(0)[i+1] - mesh->axis(0)[i])/2 + (mesh->axis(0)[1] - mesh->axis(0)[0])/2;
                            }
                        } else if (surfID == 3){
                            if (j < Ny-1){ // interior
                                jn = j+1;
                                dr = (mesh->axis(1)[j+2] - mesh->axis(1)[j])/2;
                            } else{ // periodic
                                jn = 0;
                                dr = (mesh->axis(1)[j+1] - mesh->axis(1)[j])/2 + (mesh->axis(1)[1] - mesh->axis(1)[0])/2;
                            }
                        } else{
                            if (k < Nz-1){ // interior
                                kn = k+1;
                                dr = (mesh->axis(2)[k+2] - mesh->axis(2)[k])/2;
                            } else{ // periodic
                                kn = 0;
                                dr = (mesh->axis(2)[k+1] - mesh->axis(2)[k])/2 + (mesh->axis(2)[1] - mesh->axis(2)[0])/2;
                            }
                        }

                        vars = u.matProp(in, jn, kn);
                        double uk = u(s,in,jn,kn);
                        double Tk = _converter ? uk : _converter->computeValue(vars);
                        vars[_primitive] = Tk;
                        double Dk = _mat->computeProperty(_prop, vars);
                        double dT = Dk*Tk - D*T0;

                        double flux = mesh->gradientDotN(dT/dr, i, j, k, surfID) * mesh->area(i, j, k, surfID);
                        F(s, (i*Ny + j)*Nz + k) += flux;
                        F(s, (in*Ny + jn)*Nz + kn) -= flux;
                    }

                    // Flux from boundaries
                    if (i == 0){
                        Eigen::Vector3d r{ mesh->axis(0)[i], (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(s, (i*Ny + j)*Nz + k) += u.bc(s,0)->computeFlux(u(s,i,j,k), r) * mesh->area(i,j,k,0);
                    }
                    if (i == Nx-1){
                        Eigen::Vector3d r{ mesh->axis(0)[i], (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(s, (i*Ny + j)*Nz + k) += u.bc(s,1)->computeFlux(u(s,i,j,k), r) * mesh->area(i,j,k,1);
                    }
                    if (j == 0){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, mesh->axis(1)[j], (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(s, (i*Ny + j)*Nz + k) += u.bc(s,2)->computeFlux(u(s,i,j,k), r) * mesh->area(i,j,k,2);
                    }
                    if (j == Ny-1){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, mesh->axis(1)[j], (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(s, (i*Ny + j)*Nz + k) += u.bc(s,3)->computeFlux(u(s,i,j,k), r) * mesh->area(i,j,k,3);
                    }
                    if (k == 0){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, mesh->axis(2)[k]};
                        F(s, (i*Ny + j)*Nz + k) += u.bc(s,4)->computeFlux(u(s,i,j,k), r) * mesh->area(i,j,k,4);
                    }
                    if (k == Nz-1){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, mesh->axis(2)[k]};
                        F(s, (i*Ny + j)*Nz + k) += u.bc(s,5)->computeFlux(u(s,i,j,k), r) * mesh->area(i,j,k,5);
                    }
                }
            }
        }
    }

    return F;
}

Eigen::VectorXi DiffusionKernel::stateID(const StateMesh& u) const noexcept{
    Eigen::VectorXi ind(1);
    ind(0) = u.stateID(_var);
    return ind;
}