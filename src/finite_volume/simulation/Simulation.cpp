#include "Simulation.h"
#include "Kernel.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include <iostream>

Simulation::Simulation(bool transient, std::size_t saveEveryNIterations):
    _transient(transient),
    _saveEveryNIterations(saveEveryNIterations)
{}

void Simulation::solve(StateMesh& u, double ti, double tf, double dt) const{
    _results.clear();
    auto func = [&u, this](double, const Eigen::VectorXd& arr) -> Eigen::VectorXd {
        u.flattened() = arr;
        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount());

        for (auto& kernel: _kernels){
            if (kernel){
                auto residual = kernel->computeResidual(u);
                auto indices = kernel->stateID();
                for (Eigen::Index i = 0; i < indices.size(); i++){
                    int ind = indices[i];
                    R.row(ind) += residual.row(i);
                }
            }
        }
        return Eigen::Map<Eigen::VectorXd>(R.data(), R.size());
    };

    int count = 0;
    while (true){
        if (count % _saveEveryNIterations == 0) _results[ti] = u.matrix();

        Eigen::VectorXd u0 = u.flattened();
        auto status = _integrator->integrate(func, u0, ti, dt);
        if (status != IVPStatus::Success) break;

        count++;
        if (ti >= tf) break;
        if (ti + dt > tf) ti = tf;
        else ti += dt;
    }
}