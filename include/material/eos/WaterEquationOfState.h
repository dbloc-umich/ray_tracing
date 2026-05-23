#ifndef WATER_EOS_H
#define WATER_EOS_H

#include "EquationOfState.h"

class WaterEquationOfState: public EquationOfState{
    public:
    // Mass density
    double M() const override{ return 0.018; }; // molecular mass
    double rho(double P, double T) const override; // density
    double drho_dP(double P, double T) const override; // pressure-derivative of density
    double drho_dT(double P, double T) const override; // temperature-derivative of density

    // Specific internal energy
    double E(double P, double T) const override; // specific internal energy
    double dE_dP(double P, double T) const override;
    double dE_dT(double P, double T) const override;
    double Cv(double P, double T) const override; // specific heat capacity at constant volume

    // Specific enthalpy
    double H(double P, double T) const override; // specific enthalpy
    double dH_dP(double P, double T) const override;
    double Cp(double P, double T) const override; // specific heat capacity at constant pressure

    // Transport properties
    double k(double P, double T) const override; // thermal conductivity
    double mu(double P, double T) const override; // dynamic viscosity

    // Reference states
    double Pref() const override{ return 101325.0; }; // reference pressure
    double Tref() const override{ return 273.15; }; // reference temperature

    // EOS inversion -- if P and T can be solved easily from rho and E
    bool isInvertible() const override{ return true; }; // if the inversion can be called directly from P() and T()
    double P(double rho, double E) const override; // only safe to call if isInvertible() == true
    double T(double rho, double E) const override; // only safe to call if isInvertible() == true
};

#endif