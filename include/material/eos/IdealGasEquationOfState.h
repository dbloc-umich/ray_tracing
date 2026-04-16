#ifndef IDEAL_GAS_EOS_H
#define IDEAL_GAS_EOS_H

#include "EquationOfState.h"

class IdealGasEquationOfState: public EquationOfState{
    public:
    double rho(double P, double T) const noexcept override;
    double drho_dP(double P, double T) const noexcept override{ return rho(P,T)/P; }; // pressure-derivative of density
    double drho_dT(double P, double T) const noexcept override{ return -rho(P,T)/T; }; // temperature-derivative of density
    double dCp_dT(double P, double T) const noexcept override{ return 0.0; } // temperature-derivative of heat capacity
    double H(double P, double T) const noexcept override{ return Cp(P,T)*(T-Tref()); }; // specific enthalpy
    double dH_dP(double P, double T) const noexcept override{ return 0.0; }
    double Pref() const noexcept override{ return 0.0; } // reference pressure
    double Tref() const noexcept override{ return 0.0; } // reference temperature
};

#endif