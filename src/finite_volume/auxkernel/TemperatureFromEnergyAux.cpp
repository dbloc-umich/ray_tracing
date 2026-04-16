#include "TemperatureFromEnergyAux.h"
#include "PrimitiveVariablesFromEosAux.h"

TemperatureFromEnergyAux::TemperatureFromEnergyAux(std::shared_ptr<EquationOfState> eos):
    AuxKernel(),
    _eos(eos)
{
    assert(_eos); // not null
}

double TemperatureFromEnergyAux::computeValue(const std::map<std::string, double>& u) const{
    PrimitiveVariablesFromEosAux pTSolver(_eos);
    return pTSolver.computeValue(u)[1];
}