#ifndef MATERIAL_H
#define MATERIAL_H

#include <functional>
#include <map>
#include <memory>

class EquationOfState; // for fluid properties
class MaterialProperty; // for all other properties
enum class Prop{attenuationCoefficient,
                extinctionCoefficient,
                refractiveIndex};

enum class PropVariable{concentration, pressure, temperature, velocity, wavelength};

class Material{
    public:
    using PropVars = std::map<PropVariable, double>;

    Material();
    Material(std::unique_ptr<EquationOfState> eos);
    Material(const Material&) = delete;
    Material(Material&&);
    ~Material();
    Material& operator=(const Material&) = delete;
    Material& operator=(Material&&);
    // the unique_ptr does not have access to a default deleter of an incomplete class yet,
    // so any automatically generated functions cannot be defined until the class is complete.
    
    bool hasProperty(Prop name) const noexcept;

    void addProperty(Prop name, double val) const;
    void addProperty(Prop name, std::function<double(const PropVars&)> func) const;
    void addProperty(Prop name, std::unique_ptr<MaterialProperty> prop) const;

    void replaceProperty(Prop name, double val) noexcept;
    void replaceProperty(Prop name, std::function<double(const PropVars&)> func) noexcept;
    void replaceProperty(Prop name, std::unique_ptr<MaterialProperty> prop) noexcept;

    void removeProperty(Prop name) noexcept;
    virtual double computeProperty(Prop name, const PropVars& vars = {}) const;

    // Functions from eos
    double rho(double P, double T) const noexcept; // density
    double drho_dP(double P, double T) const noexcept; // pressure-derivative of density
    double drho_dT(double P, double T) const noexcept; // temperature-derivative of density
    double beta(double P, double T) const noexcept; // thermal expansion coefficient
    double Cp(double P, double T) const noexcept; // specific heat capacity at constant pressure
    double dCp_dT(double P, double T) const noexcept; // temperature-derivative of heat capacity
    double k(double P, double T) const noexcept; // thermal conductivity
    double mu(double P, double T) const noexcept; // dynamic viscosity
    double Pr(double P, double T) const noexcept; // Prandt number
    double H(double P, double T) const noexcept; // specific enthalpy

    protected:
    std::unique_ptr<EquationOfState> _eos;
    mutable std::map<Prop, std::unique_ptr<MaterialProperty>> _props;
};

#endif