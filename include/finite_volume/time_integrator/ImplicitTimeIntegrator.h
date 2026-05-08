#ifndef IMPLICIT_TIME_INTEGRATOR_H
#define IMPLICIT_TIME_INTEGRATOR_H

#include "TimeIntegrator.h"

template<int N, int M>
class NonlinearSolver;

class ImplicitTimeIntegrator: public TimeIntegrator{
    public:
    using SolverPointer = std::unique_ptr<NonlinearSolver<Eigen::Dynamic, Eigen::Dynamic>>;
    ImplicitTimeIntegrator(SolverPointer solver = nullptr);
    ImplicitTimeIntegrator(ImplicitTimeIntegrator&&);
    virtual ~ImplicitTimeIntegrator();
    ImplicitTimeIntegrator& operator=(ImplicitTimeIntegrator&&);

    void setSolver(SolverPointer solver);

    protected:
    SolverPointer _nlSolver;   
};

#endif