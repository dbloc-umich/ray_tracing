#include "IdealGasEquationOfState.h"
#include "Constants.h"

double IdealGasEquationOfState::rho(double P, double T) const noexcept { return P*M()/(pconst::R*T); }