// #include "LaserDropletSimulation.h"
// #include "Kernel.h"
// #include "SpatialMesh.h"
// #include "StateMesh.h"
// #include "TimeIntegrator.h"

// LaserDropletSimulation::LaserDropletSimulation(bool transient, std::size_t saveEveryNIterations):
//     Simulation(transient, saveEveryNIterations)
// {}

// void LaserDropletSimulation::solve(StateMesh& u, double ti, double tf, double dt) const{
//     _results.clear();
//     auto func = [&u, this](double, const Eigen::VectorXd& arr) -> Eigen::VectorXd {
//         u.flattened() = arr;
//         Eigen::MatrixXd R = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount());

//         for (auto& kernel: _kernels){
//             if (kernel){
//                 auto residual = kernel->computeResidual(u);
//                 auto indices = kernel->stateID(u);
//                 for (Eigen::Index i = 0; i < indices.size(); i++){
//                     int ind = indices[i];
//                     R.row(ind) += residual.row(i);
//                 }
//             }
//         }
//         return Eigen::Map<Eigen::VectorXd>(R.data(), R.size());
//     };

//     auto mesh = u.mesh();
//     double m0 = 0;
//     for (Eigen::Index i = 0; i < Nr; i++){
//         for (Eigen::Index j = 0; j < Nmu; j++){
//             for (Eigen::Index k = 0; k < Nphi; k++){
//                 auto vars = u.stateMap(i,j,k);
//                 m0 += mat->computeProperty("density", vars) * mesh->volume(i,j,k);
//             }
//         }
//     }

//     int count = 0;
//     while (true){
//         if (count % _saveEveryNIterations == 0) _results[ti] = u.matrix();

//         Eigen::VectorXd u0 = u.flattened();
//         auto status = _integrator->integrate(func, u0, ti, dt);
//         if (status != IVPStatus::Success) break;

//         count++;
//         if (ti >= tf) break;
//         if (ti + dt > tf) ti = tf;
//         else ti += dt;
//     }
// }

// #include "Node.h"
// template class LaserDropletSimulation<Node>;

// #include "UNode.h"
// template class LaserDropletSimulation<UNode>;