#ifndef IDEAL_GAS_EOS_H
#define IDEAL_GAS_EOS_H

#include "EquationOfState.h"

class IdealGasEquationOfState: public EquationOfState{
    public:
    virtual ~IdealGasEquationOfState() = default;
    // Mass density
    double rho(double P, double T) const override;
    double drho_dP(double P, double T) const override; // pressure-derivative of density
    double drho_dT(double P, double T) const override; // temperature-derivative of density

    // Specific internal energy, zero pressure and constant heat capacity
    double E(double, double T) const override; // specific internal energy
    double dE_dP(double P, double T) const override{ return 0.0; }
    double dE_dT(double P, double T) const override{ return Cv(P,T); }

    // Specific enthalpy, zero pressure and constant heat capacity
    double H(double P, double T) const override; // specific enthalpy
    double dH_dP(double P, double T) const override{ return 0.0; }

    // Reference state
    double Pref() const override{ return 0.0; }
    double Tref() const override{ return 0.0; }

    // EOS inversion -- if P and T can be solved easily from rho and E
    bool isInvertible() const override{ return true; }; // if the inversion can be called directly from P() and T()
    double P(double rho, double E) const override; // only safe to call if isInvertible() == true
    double T(double rho, double E) const override; // only safe to call if isInvertible() == true
};

#endif