#ifndef MATERIAL_H
#define MATERIAL_H

#include <functional>
#include <map>
#include <memory>

class MaterialProperty; // for all other properties
enum class Prop{none,
                attenuationCoefficient,
                density,
                densityPressureDerivative,
                densityTemperatureDerivative,                
                extinctionCoefficient,
                heatCapacity,
                heatCapacityTemperatureDerivative,
                numberDensity,
                photoionizationCrossSection,
                refractiveIndex,
                thermalConductivity,
                thermalExpansionCoefficient,
                viscosity
                };

enum class PropVariable{concentration,
                        density,
                        pressure,
                        temperature,
                        velocity,
                        wavelength};

class Material{
    public:
    using PropVars = std::map<PropVariable, double>;

    Material();
    Material(const Material&) = delete;
    Material(Material&&);
    ~Material();
    Material& operator=(const Material&) = delete;
    Material& operator=(Material&&);
    // the unique_ptr does not have access to a default deleter of an incomplete class yet,
    // so any automatically generated functions cannot be defined until the class is complete.

    bool hasProperty(Prop name) const noexcept;

    void addProperty(Prop name) const;
    void addProperty(Prop name, double val) const;
    void addProperty(Prop name, std::function<double(const PropVars&)> func) const;
    void addProperty(Prop name, std::unique_ptr<MaterialProperty> prop) const;

    void replaceProperty(Prop name, double val) noexcept;
    void replaceProperty(Prop name, std::function<double(const PropVars&)> func) noexcept;
    void replaceProperty(Prop name, std::unique_ptr<MaterialProperty> prop) noexcept;

    void removeProperty(Prop name) noexcept;
    virtual double computeProperty(Prop name, const PropVars& vars = {}) const;

    // Functions from eos
    virtual double T_from_H(double H) const = 0;

    protected:
    mutable std::map<Prop, std::unique_ptr<MaterialProperty>> _props;
};

#endif