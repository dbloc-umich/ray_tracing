#include "ElectronEquationOfState.h"
#include "Constants.h"

double ElectronEquationOfState::M() const{ return pconst::m_e*pconst::N_A; }
double ElectronEquationOfState::Cv(double, double) const{ return 1.5*pconst::k_B/pconst::m_e; }
double ElectronEquationOfState::Cp(double, double) const{ return 2.5*pconst::k_B/pconst::m_e; }