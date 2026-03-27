#ifndef ELECTRON_TEMPERATURE_AUX_H
#define ELECTRON_TEMPERATURE_AUX_H

#include "AuxKernel.h"
class ElectronTemperatureAux: public AuxKernel{
    public:
    ElectronTemperatureAux(): AuxKernel() {}
    double computeValue(const std::map<std::string, double>& u) const override;
};

#endif