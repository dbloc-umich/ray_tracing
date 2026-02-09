#ifndef WATER_EOS_H
#define WATER_EOS_H

#include "EquationOfState.h"

class WaterEquationOfState: public EquationOfState{
    double rho(double P, double T) const noexcept override; // density
    double drho_dP(double P, double T) const noexcept override; // pressure-derivative of density
    double drho_dT(double P, double T) const noexcept override; // temperature-derivative of density
    double Cp(double P, double T) const noexcept override{ return 4184; } // specific heat capacity at constant pressure
    double dCp_dT(double P, double T) const noexcept override{ return 0.0; } // temperature-derivative of heat capacity
    double k(double P, double T) const noexcept override{ return 0.6065; } // thermal conductivity
    double mu(double P, double T) const noexcept override{ return 8.90e-4; } // dynamic viscosity
    double H(double P, double T) const noexcept override{ return 4184*T; } // specific enthalpy
};

#endif