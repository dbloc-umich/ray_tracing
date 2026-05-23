// The water ion of interest is H2O+, which is a neutral water molecule with one outer electron removed.

#ifndef WATER_ION_EOS_H
#define WATER_ION_EOS_H

#include "IdealGasEquationOfState.h"
class WaterIonEquationOfState: IdealGasEquationOfState{
    public:
    // Mass density
    double M() const override{ return 0.018; }; // molecular mass

    // Specific internal energy
    double Cv(double P, double T) const override; // specific heat capacity at constant volume

    // Specific enthalpy
    double Cp(double P, double T) const override; // specific heat capacity at constant pressure

    // Transport properties
    double k(double P, double T) const override{ return 0.0; } // thermal conductivity -- note: do not use k computed from this EOS. Consult WaterIonMaterialProperty.
    double mu(double P, double T) const override{ return 8.90e-4; } // dynamic viscosity
};


#endif