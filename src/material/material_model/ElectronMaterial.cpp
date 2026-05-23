#include "ElectronMaterial.h"
#include "Constants.h"
#include "ElectronEquationOfState.h"

#include <cmath>
#include <iostream>

const ElectronEquationOfState ElectronMaterial::_eos = ElectronEquationOfState();
ElectronMaterial::ElectronMaterial(const EquationOfState& nEOS):
    Material(),
    _nEOS(nEOS)
{
    addProperty("coulomb_logarithm");
    addProperty("density");
    addProperty("density_pressure_derivative");
    addProperty("density_temperature_derivative");
    addProperty("electron_ion_collision_frequency");
    addProperty("electron_neutral_collision_frequency");
    addProperty("enthalpy");
    addProperty("heat_capacity");
    addProperty("heat_capacity_temperature_derivative");
    addProperty("inverse_bremsstrahlung_frequency");
    addProperty("molecular_mass");
    addProperty("plasma_electron_frequency");
    addProperty("thermal_conductivity");
    addProperty("thermal_expansion_coefficient");
    addProperty("viscosity");
}

double ElectronMaterial::computeProperty(const std::string& name, const PropVars& vars) const{
    double ne = vars.at("electron_number_density");
    double rho = ne*pconst::m_e;
    double E = vars.at("electron_energy")/ne;
    double T = E / (1.5*pconst::k_B);
    double P = rho*(pconst::R*T)/_eos.M();
    double lambda = vars.count("wavelength") == 0 ? 589e-9 : vars.at("wavelength");
    double ni = vars.at("ion_number_density");

    if (name == "coulomb_logarithm"){
        double Z = 1.0; //computeProperty("ionization_number", vars);
        double e = pconst::e;
        double ekT = pconst::epsilon_0 * pconst::k_B * T;
        double Lambda = 4.0*mconst::pi/(Z*e*e*e) * std::sqrt(ekT*ekT*ekT/ne);
        return std::max(std::log(Lambda), 1.0);
    }
    if (name == "density") return rho;
    if (name == "density_pressure_derivative") return _eos.drho_dP(P,T);
    if (name == "density_temperature_derivative") return _eos.drho_dT(P,T);
    if (name == "electron_ion_collision_frequency"){
        double Z = 1.0; //computeProperty("ionization_number", vars);
        double lnLambda = computeProperty("coulomb_logarithm", vars);
        double e = pconst::e;
        double k = 1.0 / (4 * mconst::pi * pconst::epsilon_0);
        double kT = pconst::k_B*T;
        return (mconst::pi * ni*Z*e*e*e*e*lnLambda * k*k) / std::sqrt(pconst::m_e * (kT*kT*kT));
    }
    if (name == "electron_neutral_collision_frequency"){
        double nn = vars.at("density") / (_nEOS.M() / pconst::N_A);
        double sigma_en = 1e-19; // cross section in square meter
        double ve = std::sqrt(8*pconst::k_B*T / (mconst::pi*pconst::m_e));
        return nn*sigma_en*ve;
    }
    if (name == "enthalpy") return _eos.H(P,T);
    if (name == "heat_capacity") return _eos.Cp(P,T);
    if (name == "inverse_bremsstrahlung_frequency"){
        if (std::isnan(T)) return 0.0;
        double nu_ei = computeProperty("electron_ion_collision_frequency", vars);
        double omega = 2*mconst::pi * pconst::c/lambda;
        double omega_pe = computeProperty("plasma_electron_frequency", vars);
        return omega_pe*omega_pe / (omega*omega + nu_ei*nu_ei) * nu_ei;
    }
    if (name == "molecular_mass") return _eos.M();
    if (name == "plasma_electron_frequency") return pconst::e * std::sqrt(ne / (pconst::epsilon_0*pconst::m_e)); // angular frequency in rad/s
    if (name == "thermal_conductivity"){
        double nu = computeProperty("electron_ion_collision_frequency", vars) + computeProperty("electron_neutral_collision_frequency", vars);
        double kB = pconst::k_B;
        return 2.5*(ne*kB*kB*T)/(pconst::m_e*nu);
    }
    if (name == "thermal_expansion_coefficient") return _eos.beta(P,T);
    if (name == "viscosity") return _eos.mu(P,T);
    return Material::computeProperty(name, vars);
}

const EquationOfState& ElectronMaterial::eos() const noexcept{ return _eos; }