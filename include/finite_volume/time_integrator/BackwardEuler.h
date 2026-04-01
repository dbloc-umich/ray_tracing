#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H

#include "TimeIntegrator.h"
#include "NewtonSolver.h"

class BackwardEuler: public TimeIntegrator{
    public:
    using SolverPointer = std::unique_ptr<NonlinearSolver<Eigen::Dynamic>>;
    BackwardEuler(SolverPointer solver = std::make_unique<NewtonSolver<Eigen::Dynamic>>());   
    IVPStatus integrate(const Function& f, Eigen::VectorXd& u0, double t, double dt) const noexcept override;

    protected:
    SolverPointer _nlSolver;   
};

#endif