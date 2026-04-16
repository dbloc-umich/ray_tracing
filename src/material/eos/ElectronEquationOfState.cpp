#include "ElectronEquationOfState.h"
#include "Constants.h"

double ElectronEquationOfState::M() const noexcept{ return pconst::m_e*pconst::N_A; }
double ElectronEquationOfState::Cp(double, double) const noexcept{ return 2.5*pconst::k_B/pconst::m_e; }