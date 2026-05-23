#include "WaterIonEquationOfState.h"
#include "Constants.h"

double WaterIonEquationOfState::Cv(double, double) const{ return 3*pconst::R/M(); }
double WaterIonEquationOfState::Cp(double, double) const{ return 4*pconst::R/M(); }