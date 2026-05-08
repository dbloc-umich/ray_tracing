#include "ImplicitTimeIntegrator.h"
#include "NonlinearSolver.h"

ImplicitTimeIntegrator::ImplicitTimeIntegrator(SolverPointer solver):
    TimeIntegrator(),
    _nlSolver(std::move(solver))
{}

ImplicitTimeIntegrator::ImplicitTimeIntegrator(ImplicitTimeIntegrator&&) = default;
ImplicitTimeIntegrator::~ImplicitTimeIntegrator() = default;
ImplicitTimeIntegrator& ImplicitTimeIntegrator::operator=(ImplicitTimeIntegrator&&) = default;

void ImplicitTimeIntegrator::setSolver(SolverPointer solver){ _nlSolver = std::move(solver); }