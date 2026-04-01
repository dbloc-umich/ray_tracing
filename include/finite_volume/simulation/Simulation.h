#ifndef SIMULATION_H
#define SIMULATION_H

#include <map>
#include <memory>
#include <vector>
#include "Eigen/Dense"

class Kernel;
class StateMesh;
class TimeIntegrator;
class Simulation{
    public:
    Simulation(bool transient=true, std::size_t saveEveryNIterations=1);
    void addKernel(std::shared_ptr<Kernel> kernel) noexcept;
    void setTimeIntegrator(std::shared_ptr<TimeIntegrator> integrator) noexcept;

    virtual void solve(StateMesh& u, double ti, double tf, double dt) const;
    decltype(auto) getResults() const noexcept{ return std::move(_results); }

    protected:
    bool _transient = true; // set to true for now
    std::size_t _saveEveryNIterations;

    mutable std::map<double, Eigen::MatrixXd> _results; // full state vector results at specific times
    std::vector<std::shared_ptr<Kernel>> _kernels;
    std::shared_ptr<TimeIntegrator> _integrator;
    // std::shared_ptr<TimeStepper>

};

#endif