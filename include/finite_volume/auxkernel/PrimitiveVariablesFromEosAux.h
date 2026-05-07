// Solves for pressure and temperature using mass and energy densities information

#ifndef PRIMITIVE_VARIABLES_FROM_EOS_AUX_H
#define PRIMITIVE_VARIABLES_FROM_EOS_AUX_H

#include "AuxKernel.h"
class EquationOfState;
class PrimitiveVariablesFromEosAux: public AuxKernel<2>{
    public:
    PrimitiveVariablesFromEosAux(const EquationOfState& eos, std::string densityVar="density", std::string energyVar="energy",
                                 double pScale=1.0e-6, double tScale=1.0e-3, double rhoScale=1.0, double HScale=1.0e-3);
    Eigen::Vector2d computeValue(const std::map<std::string, double>& u) const override;

    protected:
    const EquationOfState& _eos;
    std::string _densityVar;
    std::string _energyVar;
    double _pScale; // pressure scale
    double _tScale; // temperature scale
    double _rhoScale; // density scale
    double _HScale; // enthalpy scale
};

#endif