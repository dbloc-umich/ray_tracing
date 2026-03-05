#include "Diffusion.h"
#include "MultivariateBC.h"
#include "SpatialMesh.h"
#include "StateMesh.h"
#include <iostream>

Diffusion::Diffusion(std::shared_ptr<Material> mat, PropVariable var, const Eigen::VectorXi& s, Prop prop):
    Kernel(mat, var, s),
    _prop(prop)
{}

Eigen::MatrixXd Diffusion::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::Index Nx = mesh->axisSize(0)-1; // number of cells on the 1st axis
    Eigen::Index Ny = mesh->axisSize(1)-1; // number of cells on the 2nd axis
    Eigen::Index Nz = mesh->axisSize(2)-1; // number of cells on the 3rd axis
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(1, Nx*Ny*Nz);

    if (_prop != Prop::none){
        for (Eigen::Index i = 0; i < Nx; i++){
            for (Eigen::Index j = 0; j < Ny; j++){
                for (Eigen::Index k = 0; k < Nz; k++){
                    auto vars = u.matProp(_mat, i, j, k);
                    double u0 = u(_s[0],i,j,k);
                    double D = _mat->computeProperty(_prop, vars);

                    // Flux from and to neighboring cells
                    std::vector<Eigen::Index> surfaces;
                    if (i < Nx-1 || (i == Nx-1 && u.bc(1)->isPeriodicBC())) surfaces.push_back(1);
                    if (j < Ny-1 || (j == Ny-1 && u.bc(3)->isPeriodicBC())) surfaces.push_back(3);
                    if (k < Nz-1 || (k == Nz-1 && u.bc(5)->isPeriodicBC())) surfaces.push_back(5);
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

                        double uk = u(_s[0],in,jn,kn);
                        vars = u.matProp(_mat, in, jn, kn);
                        double Dk = _mat->computeProperty(_prop, vars);
                        double du = Dk*uk - D*u0;
                        double flux = mesh->gradientDotN(du/dr, i, j, k, surfID) * mesh->area(i, j, k, surfID);
                        F(_s[0], (i*Ny + j)*Nz + k) += flux;
                        F(_s[0], (in*Ny + jn)*Nz + kn) -= flux;
                    }

                    // Flux from boundaries
                    if (i == 0){
                        Eigen::Vector3d r{ mesh->axis(0)[i], (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(_s[0], (i*Ny + j)*Nz + k) += u.bc(0)->computeFlux(u.cell(i,j,k), r)[_s[0]] * mesh->area(i,j,k,0);
                    }
                    if (i == Nx-1){
                        Eigen::Vector3d r{ mesh->axis(0)[i], (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(_s[0], (i*Ny + j)*Nz + k) += u.bc(1)->computeFlux(u.cell(i,j,k), r)[_s[0]] * mesh->area(i,j,k,1);
                    }
                    if (j == 0){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, mesh->axis(1)[j], (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(_s[0], (i*Ny + j)*Nz + k) += u.bc(2)->computeFlux(u.cell(i,j,k), r)[_s[0]] * mesh->area(i,j,k,2);
                    }
                    if (j == Ny-1){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, mesh->axis(1)[j], (mesh->axis(2)[k]+mesh->axis(2)[k+1])/2};
                        F(_s[0], (i*Ny + j)*Nz + k) += u.bc(3)->computeFlux(u.cell(i,j,k), r)[_s[0]] * mesh->area(i,j,k,3);
                    }
                    if (k == 0){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, mesh->axis(2)[k]};
                        F(_s[0], (i*Ny + j)*Nz + k) += u.bc(4)->computeFlux(u.cell(i,j,k), r)[_s[0]] * mesh->area(i,j,k,4);
                    }
                    if (k == Nz-1){
                        Eigen::Vector3d r{ (mesh->axis(0)[i]+mesh->axis(0)[i+1])/2, (mesh->axis(1)[j]+mesh->axis(1)[j+1])/2, mesh->axis(2)[k]};
                        F(_s[0], (i*Ny + j)*Nz + k) += u.bc(5)->computeFlux(u.cell(i,j,k), r)[_s[0]] * mesh->area(i,j,k,5);
                    }
                }
            }
        }
    }

    return F;
}