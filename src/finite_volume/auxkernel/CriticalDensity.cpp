// #include "CriticalDensity.h"
// #include "Constants.h"
// #include "InputParameters.h"

// CriticalDensity::CriticalDensity(const InputParameters& param):
//     AuxKernel(param),
//     _lambda(param.getRequiredScalar("laser_wavelength"))
// {}

// double CriticalDensity::computeValue(const ConstCell& u) const{
//     double m_e = pconst::m_e;
//     double c = pconst::c;
//     double e = pconst::e;
//     double pi = mconst::pi;
//     return m_e*pi*c*c/(_lambda*_lambda*e*e);
// }

