#include "WaterMaterial.h"
#include "Constants.h"
#include "ElectronTemperatureAux.h"
#include "TemperatureFromEnergyAux.h"
#include "WaterEquationOfState.h"
#include <iostream>

const std::shared_ptr<WaterEquationOfState> WaterMaterial::_eos = std::make_shared<WaterEquationOfState>();
WaterMaterial::WaterMaterial():
    Material()
{
    addProperty("2-photon_ionization_rate_coefficient");
    addProperty("3-photon_ionization_rate_coefficient");
    addProperty("4-photon_ionization_rate_coefficient");
    addProperty("attenuation_coefficient");
    addProperty("binding_energy");
    addProperty("coulomb_logarithm");
    addProperty("density");
    addProperty("density_pressure_derivative");
    addProperty("density_temperature_derivative");
    addProperty("electron_density");
    addProperty("electron_ion_collision_frequency");
    addProperty("electron_neutral_collision_frequency");
    addProperty("enthalpy");
    addProperty("extinction_coefficient");
    addProperty("heat_capacity");
    addProperty("heat_capacity_temperature_derivative");
    addProperty("inverse_bremsstrahlung_frequency");
    addProperty("ionization_number");
    addProperty("molecular_mass");
    addProperty("number_density");
    addProperty("quantum_yield");
    addProperty("refractive_index");
    addProperty("thermal_conductivity");
    addProperty("thermal_expansion_coefficient");
    addProperty("viscosity");
}

double WaterMaterial::computeProperty(const std::string& name, const PropVars& vars) const{
    double P = vars.count("pressure") == 0 ? _eos->Pref() : vars.at("pressure");
    double lambda = vars.count("wavelength") == 0 ? 589e-9 : vars.at("wavelength");
    double rho = vars.count("density") == 0 ? _eos->rho(_eos->Pref(), _eos->Tref()) : vars.at("density");
    double /*H,*/ T;
    if (vars.count("energy") == 0){
        if (vars.count("temperature") == 0){
            T = _eos->Tref();
            // H = 0.0;
        } else{
            T = vars.at("temperature");
            // H = _eos->H(P, T);
        }
    } else{
        // H = vars.at("energy") / rho;
        TemperatureFromEnergyAux aux(_eos);
        T = aux.computeValue(vars);
    }
    double ni = vars.count("ion_number_density") == 0 ? 0 : vars.at("ion_number_density");
    // double Ee = vars.count("electron_energy") == 0 ? 0.0 : vars.at("electron_energy") / ni;
    double Te = ElectronTemperatureAux().computeValue(vars);

    if (name == "2-photon_ionization_rate_coefficient") return 1.0e-52 * 1.0e-8;  // lit data = 1.0e-52 cm4 s, convert to m4 s
    if (name == "3-photon_ionization_rate_coefficient") return 1.0e-84 * 1.0e-12;  // lit data = 1.0e-84 cm6 s2, convert to m6 s2
    if (name == "4-photon_ionization_rate_coefficient") return 1.0e-118 * 1.0e-16; // lit data = 1.0e-118 cm8 s3, convert to m8 s3
    if (name == "attenuation_coefficient"){
        double kappa = computeProperty("extinction_coefficient", vars);
        return 4*mconst::pi*kappa/lambda;
    }
    if (name == "binding_energy") return 5.9280535458e-19; // 3.7 eV
    if (name == "coulomb_logarithm"){
        if (std::isnan(Te)) return 0.0;
        double Z = computeProperty("ionization_number", vars);
        double e = pconst::e;
        double kT = pconst::k_B * Te;
        double ne = computeProperty("electron_density", vars);
        double Lambda = 3.0/(2*Z*e*e*e) * std::sqrt(kT*kT*kT / (mconst::pi*ne));
        return std::log(Lambda);
    }
    if (name == "density") return rho;
    if (name == "density_pressure_derivative") return _eos->drho_dP(P,T);
    if (name == "density_temperature_derivative") return _eos->drho_dT(P,T);
    if (name == "electron_density") return ni * computeProperty("ionization_number", vars);
    if (name == "electron_ion_collision_frequency"){
        if (std::isnan(Te)) return 0.0;
        double Z = computeProperty("ionization_number", vars);
        double lnLambda = computeProperty("coulomb_logarithm", vars);
        double ne = computeProperty("electron_density", vars);
        double e = pconst::e;
        return 4.0/3 * std::sqrt(2*mconst::pi/pconst::m_e) * ne*Z*e*e*e*e*lnLambda * std::pow(pconst::k_B*Te, -1.5);
    }
    if (name == "electron_neutral_collision_frequency"){
        if (std::isnan(Te)) return 0.0;
        double nn = computeProperty("number_density", vars) - ni;
        double sigma_en = 1e-19; // square meter
        double ve = std::sqrt(8*pconst::k_B*Te / (mconst::pi*pconst::m_e));
        return nn*sigma_en*ve;
    }
    if (name == "enthalpy") return _eos->H(P,T);
    if (name == "extinction_coefficient"){
        double rho_bar = _eos->rho(P,T) / _eos->rho(101325, 293.15); 
        if (lambda <= _lambda[0]) return _kappa[0]*rho_bar;
        if (lambda >= _lambda[_lambda.size()-1]) return _kappa[_kappa.size()-1]*rho_bar;
        std::size_t ind = std::upper_bound(_lambda.cbegin(), _lambda.cend(), lambda) - _lambda.cbegin();
        return (_kappa[ind] - (_kappa[ind]-_kappa[ind-1])/(_lambda[ind]-_lambda[ind-1]) * (_lambda[ind]-lambda))*rho_bar;
    }
    if (name == "heat_capacity") return _eos->Cp(P,T);
    if (name == "heat_capacity_temperature_derivative") return _eos->dCp_dT(P,T);
    if (name == "inverse_bremsstrahlung_frequency"){
        if (std::isnan(Te)) return 0.0;
        double nu_ei = computeProperty("electron_ion_collision_frequency", vars);
        double v = pconst::c / computeProperty("refractive_index", vars);
        double omega = 2*mconst::pi * v/lambda;
        return omega*omega / (omega*omega + nu_ei*nu_ei) * nu_ei;
    }
    if (name == "ionization_number") return 1.0;
    if (name == "molecular_mass") return _eos->M();
    if (name == "number_density") return _eos->rho(P,T) / _eos->M() * pconst::N_A;
    if (name == "photoionizationCrossSection") return 6.30e+6; // in barns
    if (name == "quantum_yield") return 0.9;
    if (name == "refractive_index"){
        double T_bar = T/273.15;
        double rho_bar = _eos->rho(P,T) / 1000;
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

    if (name == "thermal_conductivity") return _eos->k(P,T);
    if (name == "thermal_expansion_coefficient") return _eos->beta(P,T);
    if (name == "viscosity") return _eos->mu(P,T);
    return Material::computeProperty(name, vars);
}

double WaterMaterial::T_from_H(double H) const{ return _eos->T(H); }