#include "WaterIonEquationOfState.h"
#include "Constants.h"

double WaterIonEquationOfState::Cp(double, double) const noexcept{ return 1.4*pconst::R/M(); }