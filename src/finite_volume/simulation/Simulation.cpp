#include "Simulation.h"
#include "Kernel.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include <iostream>

Simulation::Simulation(bool transient, std::size_t saveEveryNIterations):
    _transient(transient),
    _saveEveryNIterations(saveEveryNIterations)
{}

void Simulation::addKernel(std::shared_ptr<Kernel> kernel) noexcept{
    if (kernel) _kernels.push_back(kernel);
}

void Simulation::setTimeIntegrator(std::shared_ptr<TimeIntegrator> integrator) noexcept{
    // Checks against setting a null integrator to a transient simulation
    if (_transient || integrator) _integrator = integrator;
}

void Simulation::solve(StateMesh& u, double ti, double tf, double dt) const{
    _results.clear();
    auto func = [&u, this](double, const Eigen::VectorXd& arr) -> Eigen::VectorXd {
        u.flattened() = std::move(arr);
        // std::cout << "Initial matrix: " << std::endl << u.matrix() << std::endl;
        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount());

        for (auto& kernel: _kernels){
            if (kernel){
                auto residual = kernel->computeResidual(u);
                auto indices = kernel->stateID(u);
                for (Eigen::Index i = 0; i < indices.size(); i++){
                    int ind = indices[i];
                    R.row(ind) += residual.row(i);
                }
            }
        }
        // std::cout << "Residual: " << std::endl << R << std::endl << std::endl;
        return Eigen::Map<Eigen::VectorXd>(R.data(), R.size());
    };

    Eigen::VectorXd u0 = u.flattened();
    int count = 0;
    if (_integrator){
        while (true){
            if (count % _saveEveryNIterations == 0) _results[ti] = u.matrix();
            double step = std::min(tf-ti, dt);
            
            auto status = _integrator->integrate(func, u0, ti, step);
            if (status != IVPStatus::Success){
                std::cerr << "ERROR: Time integration failed." << std::endl;
                break;
            }

            count++;
            if (ti >= tf) break;
            ti += step;
        }
    }
}