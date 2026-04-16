// Solves for pressure and temperature using mass and energy densities information

#ifndef PRIMITIVE_VARIABLES_FROM_EOS_AUX_H
#define PRIMITIVE_VARIABLES_FROM_EOS_AUX_H

#include "AuxKernel.h"
class EquationOfState;
class PrimitiveVariablesFromEosAux: public AuxKernel<2>{
    public:
    PrimitiveVariablesFromEosAux(std::shared_ptr<EquationOfState> eos);
    Eigen::Vector2d computeValue(const std::map<std::string, double>& u) const override;

    protected:
    std::shared_ptr<EquationOfState> _eos;
};

#endif