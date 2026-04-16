#include "PrimitiveVariablesFromEosAux.h"
#include "EquationOfState.h"
#include "NewtonSolver.h"
#include <algorithm>

PrimitiveVariablesFromEosAux::PrimitiveVariablesFromEosAux(std::shared_ptr<EquationOfState> eos):
    AuxKernel(),
    _eos(eos)
{
    assert(_eos); // not null
}

Eigen::Vector2d PrimitiveVariablesFromEosAux::computeValue(const std::map<std::string, double>& u) const{
    auto func = [this, &u](const Eigen::Vector2d& pT) -> Eigen::Vector2d {
        double p = pT[0];
        double T = pT[1];
        double rho = _eos->rho(p,T);
        double H = _eos->H(p, T);

        Eigen::Vector2d f;
        f[0] = rho - u.at("density"); // density equation: rho = rho(p,T)
        f[1] = rho*H - (u.at("energy") + p); // energy equation: rho*E + p = rho(p,T)*H(p,T)
        return f;
    };

    auto jac = [this, &u](const Eigen::Vector2d& pT) -> Eigen::Matrix2d {
        double p = pT[0];
        double T = pT[1];
        double rho = _eos->rho(p,T);
        double drho_dP = _eos->drho_dP(p,T);
        double drho_dT = _eos->drho_dT(p,T);
        double H = _eos->H(p, T);
        double Cp = _eos->Cp(p,T);
        double dH_dP = _eos->dH_dP(p,T);

        Eigen::Matrix2d J;
        J(0,0) = drho_dP;
        J(0,1) = drho_dT;
        J(1,0) = drho_dP*H + rho*dH_dP - 1.0;
        J(1,1) = drho_dT*H + rho*Cp;
        return J;
    };

    Eigen::Vector2d pT{_eos->Pref(), _eos->Tref()};
    NewtonSolver<2> newton(func, jac);
    auto status = newton.solve(pT);
    switch (status){
        case (NLStatus::Success):
            return pT;
        default:
            throw std::runtime_error("ERROR: In PrimitiveVariablesFromEosAux - Unable to solve for pressure and temperature.");
    }
}