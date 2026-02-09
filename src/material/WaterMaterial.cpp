#include "WaterMaterial.h"
#include "Constants.h"
#include "WaterEquationOfState.h"

WaterMaterial::WaterMaterial():
    Material(std::make_unique<WaterEquationOfState>())
{
    addProperty(Prop::refractiveIndex, 0);
    addProperty(Prop::extinctionCoefficient, 0);
    addProperty(Prop::attenuationCoefficient, 0);
}

double WaterMaterial::computeProperty(Prop name, const PropVars& vars) const{
    double P = vars.count(PropVariable::pressure) == 0 ? 101325 : vars.at(PropVariable::pressure);
    double T = vars.count(PropVariable::temperature) == 0 ? 300 : vars.at(PropVariable::temperature);
    double lambda = vars.count(PropVariable::wavelength) == 0 ? 589e-9 : vars.at(PropVariable::wavelength);

    switch(name){
        case Prop::refractiveIndex: {
            double T_bar = T/273.15;
            double rho_bar = rho(P,T) / 1000;
            double lam_bar = lambda/589e-9;
            double rhs = 0.244257733;
            rhs += 0.00974634476*rho_bar;
            rhs -= 0.00373234996*T_bar;
            rhs += 0.000268678472*lam_bar*lam_bar*T_bar;
            rhs += 0.0015892057/(lam_bar*lam_bar);
            rhs += 0.00245934259/(lam_bar*lam_bar - 0.229202*0.229202);
            rhs += 0.90070492/(lam_bar*lam_bar - 5.432937*5.432937);
            rhs -= 0.0166626219*rho_bar*rho_bar;
            rhs *= rho_bar;
            return std::sqrt((2*rhs+1)/(1-rhs));
        }
        case Prop::extinctionCoefficient: {
            double rho_bar = rho(P,T) / rho(101325, 293.15); 
            if (lambda <= _lambda[0]) return _kappa[0]*rho_bar;
            if (lambda >= _lambda[_lambda.size()-1]) return _kappa[_kappa.size()-1]*rho_bar;
            std::size_t ind = std::upper_bound(_lambda.cbegin(), _lambda.cend(), lambda) - _lambda.cbegin();
            return (_kappa[ind] - (_kappa[ind]-_kappa[ind-1])/(_lambda[ind]-_lambda[ind-1]) * (_lambda[ind]-lambda))*rho_bar;
        }
        case Prop::attenuationCoefficient: {
            double kappa = computeProperty(Prop::extinctionCoefficient, vars);
            return 4*mconst::pi*kappa/lambda;
        }
        default:
            return Material::computeProperty(name, vars);
    }
}