// This class's implementation is derived from PrimitiveVariablesFromEosAux, where the pressure is discarded and only the temperature is needed.

#ifndef TEMPERATURE_FROM_ENERGY_AUX_H
#define TEMPERATURE_FROM_ENERGY_AUX_H

#include "AuxKernel.h"
#include "PrimitiveVariablesFromEosAux.h"
class EquationOfState;
class TemperatureFromEnergyAux: public AuxKernel<1>{
    public:
    TemperatureFromEnergyAux(const EquationOfState& eos, std::string densityVar="density", std::string energyVar="energy",
                             double pScale=1.0e-6, double tScale=1.0e-3, double rhoScale=1.0, double HScale=1.0e-3);
    double computeValue(const std::map<std::string, double>& u) const override;

    protected:
    PrimitiveVariablesFromEosAux _ptAux;
};

#endif