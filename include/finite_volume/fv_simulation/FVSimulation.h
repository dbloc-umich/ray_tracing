#ifndef FV_SIMULATION_H
#define FV_SIMULATION_H

#include <memory>
#include <vector>

class FVKernel;
class FVStateMesh;
class FVTimeIntegrator;
class FVSimulation{
    protected:
    std::vector<std::shared_ptr<FVKernel>> _kernels;
    std::shared_ptr<FVTimeIntegrator> _integrator;
};

#endif