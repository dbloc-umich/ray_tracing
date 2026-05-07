#include "IdealGasEquationOfState.h"
#include "Constants.h"

double IdealGasEquationOfState::rho(double P, double T) const noexcept { return P*M()/(pconst::R*T); }
double IdealGasEquationOfState::drho_dP(double P, double T) const noexcept { return M()/(pconst::R*T); }
double IdealGasEquationOfState::drho_dT(double P, double T) const noexcept { return -P*M()/(pconst::R*T*T); }