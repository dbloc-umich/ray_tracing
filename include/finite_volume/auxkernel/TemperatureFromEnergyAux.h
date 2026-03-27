#ifndef TEMPERATURE_FROM_ENERGY_AUX_H
#define TEMPERATURE_FROM_ENERGY_AUX_H

#include "AuxKernel.h"
class EquationOfState;
class TemperatureFromEnergyAux: public AuxKernel{
    public:
    TemperatureFromEnergyAux(std::shared_ptr<EquationOfState> eos);
    double computeValue(const std::map<std::string, double>& u) const override;

    protected:
    std::shared_ptr<EquationOfState> _eos;
};

#endif