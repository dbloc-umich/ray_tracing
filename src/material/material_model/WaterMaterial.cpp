#include "WaterMaterial.h"
#include "Constants.h"
#include "WaterEquationOfState.h"
#include <iostream>

const WaterEquationOfState WaterMaterial::_eos = WaterEquationOfState();
WaterMaterial::WaterMaterial():
    Material()
{
    addProperty("2-photon_ionization_rate_coefficient");
    addProperty("3-photon_ionization_rate_coefficient");
    addProperty("4-photon_ionization_rate_coefficient");
    addProperty("attenuation_coefficient");
    addProperty("binding_energy");
    addProperty("density");
    addProperty("density_pressure_derivative");
    addProperty("density_temperature_derivative");
    addProperty("enthalpy");
    addProperty("extinction_coefficient");
    addProperty("heat_capacity");
    addProperty("molecular_mass");
    addProperty("number_density");
    addProperty("quantum_yield");
    addProperty("refractive_index");
    addProperty("thermal_conductivity");
    addProperty("thermal_expansion_coefficient");
    addProperty("viscosity");
}

double WaterMaterial::computeProperty(const std::string& name, const PropVars& vars) const{
    double lambda = vars.count("wavelength") == 0 ? 589e-9 : vars.at("wavelength");
    double rho, E, P, T;
    if (vars.count("density") == 0 || vars.count("energy") == 0){
        P = _eos.Pref();
        T = _eos.Tref();
        rho = _eos.rho(P, T);
        E = _eos.E(P, T);
    } else{
        rho = vars.at("density");
        E = vars.at("energy")/rho;
        P = _eos.P(rho, E);
        T = _eos.T(rho, E);
    }

    if (name == "2-photon_ionization_rate_coefficient") return 1.0e-52 * 1.0e-8;  // lit data = 1.0e-52 cm4 s, convert to m4 s
    if (name == "3-photon_ionization_rate_coefficient") return 1.0e-84 * 1.0e-12;  // lit data = 1.0e-84 cm6 s2, convert to m6 s2
    if (name == "4-photon_ionization_rate_coefficient") return 1.0e-118 * 1.0e-16; // lit data = 1.0e-118 cm8 s3, convert to m8 s3
    if (name == "attenuation_coefficient"){
        double kappa = computeProperty("extinction_coefficient", vars);
        return 4*mconst::pi*kappa/lambda;
    }
    if (name == "binding_energy") return 5.9280535458e-19; // 3.7 eV
    if (name == "density") return rho;
    if (name == "density_pressure_derivative") return _eos.drho_dP(P,T);
    if (name == "density_temperature_derivative") return _eos.drho_dT(P,T);
    if (name == "enthalpy") return _eos.H(P,T);
    if (name == "extinction_coefficient"){
        double rho_bar = rho / 1000; 
        if (lambda <= _lambda[0]) return _kappa[0]*rho_bar;
        if (lambda >= _lambda[_lambda.size()-1]) return _kappa[_kappa.size()-1]*rho_bar;
        std::size_t ind = std::upper_bound(_lambda.cbegin(), _lambda.cend(), lambda) - _lambda.cbegin();
        return (_kappa[ind] - (_kappa[ind]-_kappa[ind-1])/(_lambda[ind]-_lambda[ind-1]) * (_lambda[ind]-lambda))*rho_bar;
    }
    if (name == "heat_capacity") return _eos.Cp(P,T);
    if (name == "molecular_mass") return _eos.M();
    if (name == "number_density") return _eos.rho(P,T) / _eos.M() * pconst::N_A;
    if (name == "orbital_kinetic_energy") return pconst::Ry; // assumed to just be hydrogen for now
    if (name == "photoionization_cross_section") return 6.30e+6; // in barns
    if (name == "quantum_yield") return 0.9;
    if (name == "refractive_index"){
        double T_bar = T/273.15;
        double rho_bar = rho / 1000;
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

    if (name == "thermal_conductivity") return _eos.k(P,T);
    if (name == "thermal_expansion_coefficient") return _eos.beta(P,T);
    if (name == "viscosity") return _eos.mu(P,T);
    return Material::computeProperty(name, vars);
}

const EquationOfState& WaterMaterial::eos() const noexcept{ return _eos; }