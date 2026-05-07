#include "PrimitiveVariablesFromEosAux.h"
#include "EquationOfState.h"
#include "NewtonSolver.h"
#include <algorithm>

PrimitiveVariablesFromEosAux::PrimitiveVariablesFromEosAux(const EquationOfState& eos, std::string densityVar, std::string energyVar,
                                                           double pScale, double tScale, double rhoScale, double HScale):
    AuxKernel(),
    _eos(eos),
    _densityVar(std::move(densityVar)),
    _energyVar(std::move(energyVar)),
    _pScale(pScale),
    _tScale(tScale),
    _rhoScale(rhoScale),
    _HScale(HScale)
{
    assert(pScale > 0.0);
    assert(tScale > 0.0);
    assert(rhoScale > 0.0);
    assert(HScale > 0.0);
}

Eigen::Vector2d PrimitiveVariablesFromEosAux::computeValue(const std::map<std::string, double>& u) const{
    auto func = [this, &u](const Eigen::Vector2d& PThat) -> Eigen::Vector2d {
        double P = PThat[0] / _pScale;
        double T = PThat[1] / _tScale;
        double rho = _eos.rho(P,T);
        double H = _eos.H(P, T);

        Eigen::Vector2d f;
        f[0] = (rho - u.at(_densityVar)) * _rhoScale; // density equation: rho = rho(p,T)
        f[1] = (rho*H - (u.at(_energyVar) + P)) * _rhoScale*_HScale; // energy equation: rho*E + p = rho(p,T)*H(p,T)
        return f;
    };

    auto jac = [this, &u](const Eigen::Vector2d& PThat) -> Eigen::Matrix2d {
        double P = PThat[0] / _pScale;
        double T = PThat[1] / _tScale;
        double rho = _eos.rho(P,T);
        double drho_dPhat = _eos.drho_dP(P,T) / _pScale;
        double drho_dThat = _eos.drho_dT(P,T) / _tScale;
        double H = _eos.H(P, T);
        double dH_dThat = _eos.Cp(P,T) / _tScale;
        double dH_dPhat = _eos.dH_dP(P,T) / _pScale;

        Eigen::Matrix2d J;
        J(0,0) = drho_dPhat * _rhoScale;
        J(0,1) = drho_dThat * _rhoScale;
        J(1,0) = (drho_dPhat*H + rho*dH_dPhat - 1.0/_pScale) * _rhoScale*_HScale;
        J(1,1) = (drho_dThat*H + rho*dH_dThat) * _rhoScale*_HScale;
        return J;
    };

    Eigen::Vector2d PThat{_eos.Pref()*_pScale, _eos.Tref()*_tScale};
    NewtonSolver<2> newton(func, jac);
    auto status = newton.solve(PThat);
    switch (status){
        case (NLStatus::Success):
            PThat[0] /= _pScale;
            PThat[1] /= _tScale;
            return PThat;
        case (NLStatus::SingularityError):
            throw std::runtime_error("ERROR: In PrimitiveVariablesFromEosAux - Singular Jacobian found.");
        default:
            throw std::runtime_error("ERROR: In PrimitiveVariablesFromEosAux - Unable to solve for pressure and temperature.");
    }
}