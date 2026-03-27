#include "TemperatureFromEnergyAux.h"
#include "EquationOfState.h"
#include <algorithm>

TemperatureFromEnergyAux::TemperatureFromEnergyAux(std::shared_ptr<EquationOfState> eos):
    AuxKernel(),
    _eos(eos)
{
    assert(_eos); // not null
}

double TemperatureFromEnergyAux::computeValue(const std::map<std::string, double>& u) const{
    if (u.count("energy") == 0) throw std::out_of_range("ERROR: Required variable(s) not found.");
    double rhoH = u.at("energy");
    double rho = (u.count("density") == 0) ? _eos->rho(_eos->Pref(), _eos->Tref()) : u.at("density");
    return _eos->T(rhoH/rho);
}