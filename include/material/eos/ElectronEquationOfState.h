#ifndef ELECTRON_EOS_H
#define ELECTRON_EOS_H

#include "IdealGasEquationOfState.h"

class ElectronEquationOfState: public IdealGasEquationOfState{
    public:
    double M() const noexcept override; // molecular mass
    double Cp(double P, double T) const noexcept override; // specific heat capacity at constant pressure
    double k(double P, double T) const noexcept override{ return 0.0; } // thermal conductivity -- note: do not use k computed from this EOS. Consult ElectronMaterialProperty.
    double mu(double P, double T) const noexcept override{ return 0.0; } // dynamic viscosity
    double Pref() const noexcept override{ return 1.0; } // reference pressure
    double Tref() const noexcept override{ return 1.0; } // reference temperature
};
#endif