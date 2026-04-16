// The water ion of interest is H2O+, which is a neutral water molecule with one outer electron removed.

#ifndef WATER_ION_EOS_H
#define WATER_ION_EOS_H

#include "IdealGasEquationOfState.h"
class WaterIonEquationOfState: IdealGasEquationOfState{
    public:
    double M() const noexcept override{ return 0.018; }; // molecular mass
    double Cp(double P, double T) const noexcept override; // specific heat capacity at constant pressure
    double k(double P, double T) const noexcept override{ return 0.0; } // thermal conductivity -- note: do not use k computed from this EOS. Consult WaterIonMaterialProperty.
    double mu(double P, double T) const noexcept override{ return 8.90e-4; } // dynamic viscosity
};


#endif