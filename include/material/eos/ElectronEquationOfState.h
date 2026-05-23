#ifndef ELECTRON_EOS_H
#define ELECTRON_EOS_H

#include "IdealGasEquationOfState.h"

class ElectronEquationOfState: public IdealGasEquationOfState{
    public:
    // Mass density
    double M() const override; // molecular mass

    // Specific internal energy
    double Cv(double P, double T) const override; // specific heat capacity at constant volume

    // Specific enthalpy
    double Cp(double P, double T) const override; // specific heat capacity at constant pressure

    // Transport properties
    double k(double P, double T) const override{ return 0.0; }; // thermal conductivity -- note: do not use k computed from this EOS. Consult ElectronMaterialProperty.
    double mu(double P, double T) const override{ return 0.0; }; // dynamic viscosity
};
#endif