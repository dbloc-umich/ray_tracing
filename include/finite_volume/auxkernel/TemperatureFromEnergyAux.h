// This class's implementation is derived from PrimitiveVariablesFromEosAux, where the pressure is discarded and only the temperature is needed.

#ifndef TEMPERATURE_FROM_ENERGY_AUX_H
#define TEMPERATURE_FROM_ENERGY_AUX_H

#include "AuxKernel.h"
class EquationOfState;
class TemperatureFromEnergyAux: public AuxKernel<1>{
    public:
    TemperatureFromEnergyAux(std::shared_ptr<EquationOfState> eos);
    double computeValue(const std::map<std::string, double>& u) const override;

    protected:
    std::shared_ptr<EquationOfState> _eos;
};

#endif