#ifndef EQUATION_OF_STATE_H
#define EQUATION_OF_STATE_H

class EquationOfState{
    public:
    virtual ~EquationOfState() = default;

    // Mass density
    virtual double M() const = 0; // molecular mass
    virtual double rho(double P, double T) const = 0; // density
    virtual double drho_dP(double P, double T) const = 0; // pressure-derivative of density
    virtual double drho_dT(double P, double T) const = 0; // temperature-derivative of density
    virtual double beta(double P, double T) const{ return -drho_dT(P,T)/rho(P,T); } // thermal expansion coefficient

    // Specific internal energy
    virtual double E(double P, double T) const = 0; // specific internal energy
    virtual double dE_dP(double P, double T) const = 0;
    virtual double dE_dT(double P, double T) const = 0;
    virtual double Cv(double P, double T) const = 0; // specific heat capacity at constant volume

    // Specific enthalpy
    virtual double H(double P, double T) const = 0; // specific enthalpy
    virtual double dH_dP(double P, double T) const = 0;
    virtual double dH_dT(double P, double T) const{ return Cp(P,T); }
    virtual double Cp(double P, double T) const = 0; // specific heat capacity at constant pressure

    // Transport properties
    virtual double k(double P, double T) const = 0; // thermal conductivity
    virtual double mu(double P, double T) const = 0; // dynamic viscosity
    virtual double Pr(double P, double T) const noexcept{ return mu(P,T)*Cp(P,T)/k(P,T); } // Prandt number

    // Reference states
    virtual double Pref() const = 0; // reference pressure
    virtual double Tref() const = 0; // reference temperature

    // EOS inversion -- if P and T can be solved easily from rho and E
    virtual bool isInvertible() const = 0; // if the inversion can be called directly from P() and T()
    virtual double P(double rho, double E) const = 0; // only safe to call if isInvertible() == true
    virtual double T(double rho, double E) const = 0; // only safe to call if isInvertible() == true
};

#endif