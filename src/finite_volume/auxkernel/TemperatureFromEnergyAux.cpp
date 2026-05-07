#include "TemperatureFromEnergyAux.h"

TemperatureFromEnergyAux::TemperatureFromEnergyAux(const EquationOfState& eos, std::string densityVar, std::string energyVar,
                                                   double pScale, double tScale, double rhoScale, double HScale):
    AuxKernel(),
    _ptAux(eos, densityVar, energyVar, pScale, tScale, rhoScale, HScale)
{}

double TemperatureFromEnergyAux::computeValue(const std::map<std::string, double>& u) const{ return _ptAux.computeValue(u)[1]; }