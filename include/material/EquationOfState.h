#ifndef EQUATION_OF_STATE_H
#define EQUATION_OF_STATE_H

class EquationOfState{
    public:
    virtual ~EquationOfState() = default;

    // Required fluid properties
    virtual double rho(double P, double T) const noexcept = 0; // density
    virtual double drho_dP(double P, double T) const noexcept = 0; // pressure-derivative of density
    virtual double drho_dT(double P, double T) const noexcept = 0; // temperature-derivative of density
    virtual double beta(double P, double T) const noexcept{ return -drho_dT(P,T)/rho(P,T); } // thermal expansion coefficient
    virtual double Cp(double P, double T) const noexcept = 0; // specific heat capacity at constant pressure
    virtual double dCp_dT(double P, double T) const noexcept = 0; // temperature-derivative of heat capacity
    virtual double k(double P, double T) const noexcept = 0; // thermal conductivity
    virtual double mu(double P, double T) const noexcept = 0; // dynamic viscosity
    virtual double Pr(double P, double T) const noexcept{ return Cp(P,T)/(mu(P,T)*k(P,T)); } // Prandt number
    virtual double H(double P, double T) const noexcept = 0; // specific enthalpy
};

#endif