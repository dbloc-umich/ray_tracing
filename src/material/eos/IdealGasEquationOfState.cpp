#include "IdealGasEquationOfState.h"
#include "Constants.h"

double IdealGasEquationOfState::rho(double P, double T) const{ return P*M()/(pconst::R*T); }
double IdealGasEquationOfState::drho_dP(double P, double T) const{ return M()/(pconst::R*T); }
double IdealGasEquationOfState::drho_dT(double P, double T) const{ return -P*M()/(pconst::R*T*T); }

double IdealGasEquationOfState::E(double P, double T) const{ return Cv(P,T) * (T-Tref()); }
double IdealGasEquationOfState::H(double P, double T) const{ return pconst::R*Tref() + Cp(P,T) * (T-Tref()); }

double IdealGasEquationOfState::P(double rho, double E) const{ return rho*(pconst::R*T(rho, E)) / M(); }
double IdealGasEquationOfState::T(double rho, double E) const{ return E/Cv(0,0) + Tref(); }