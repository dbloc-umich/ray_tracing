#include "WaterEquationOfState.h"
#include "CoolProp.h"

double WaterEquationOfState::rho(double P, double T) const{ return CoolProp::PropsSI("Dmass", "P", P, "T", T, "Water"); }
double WaterEquationOfState::drho_dP(double P, double T) const{ return CoolProp::PropsSI("d(Dmass)/d(P)|T", "P", P, "T", T, "Water"); }
double WaterEquationOfState::drho_dT(double P, double T) const{ return CoolProp::PropsSI("d(Dmass)/d(T)|P", "P", P, "T", T, "Water"); }

double WaterEquationOfState::E(double P, double T) const{ return CoolProp::PropsSI("Umass", "P", P, "T", T, "Water"); }
double WaterEquationOfState::dE_dP(double P, double T) const{ return CoolProp::PropsSI("d(Umass)/d(P)|T", "P", P, "T", T, "Water"); }
double WaterEquationOfState::dE_dT(double P, double T) const{ return CoolProp::PropsSI("d(Umass)/d(T)|P", "P", P, "T", T, "Water"); }
double WaterEquationOfState::Cv(double P, double T) const{ return CoolProp::PropsSI("Cvmass", "P", P, "T", T, "Water"); }

double WaterEquationOfState::H(double P, double T) const{ return CoolProp::PropsSI("Hmass", "P", P, "T", T, "Water"); }
double WaterEquationOfState::dH_dP(double P, double T) const{ return CoolProp::PropsSI("d(Hmass)/d(P)|T", "P", P, "T", T, "Water"); }
double WaterEquationOfState::Cp(double P, double T) const{ return CoolProp::PropsSI("d(Hmass)/d(T)|P", "P", P, "T", T, "Water"); }

double WaterEquationOfState::k(double P, double T) const{ return CoolProp::PropsSI("conductivity", "P", P, "T", T, "Water"); }
double WaterEquationOfState::mu(double P, double T) const{ return CoolProp::PropsSI("viscosity", "P", P, "T", T, "Water"); }

double WaterEquationOfState::P(double rho, double E) const{ return CoolProp::PropsSI("P", "Dmass", rho, "Umass", E, "Water"); }
double WaterEquationOfState::T(double rho, double E) const{ return CoolProp::PropsSI("T", "Dmass", rho, "Umass", E, "Water"); }