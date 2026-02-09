#include "Material.h"
#include "ConstantProperty.h"
#include "EquationOfState.h"
#include "FunctionProperty.h"

#include <algorithm>
#include <stdexcept>

Material::Material():
    _eos(nullptr)
{}

Material::Material(std::unique_ptr<EquationOfState> eos):
    _eos(std::move(eos))
{}

Material::~Material() = default;
Material::Material(Material&&) = default;
Material& Material::operator=(Material&&) = default;

bool Material::hasProperty(Prop name) const noexcept{ return _props.count(name); }

void Material::addProperty(Prop name, double val) const{
    if (hasProperty(name))
        throw std::runtime_error("ERROR: A property of the same name already exists.");
    _props[name] = std::make_unique<ConstantProperty>(val);
    _props[name]->setMaterial(this);
}

void Material::addProperty(Prop name, std::function<double(const PropVars&)> func) const{
    if (hasProperty(name))
        throw std::runtime_error("ERROR: A property of the same name already exists.");
    _props[name] = std::make_unique<FunctionProperty>(func);
    _props[name]->setMaterial(this);
}

void Material::addProperty(Prop name, std::unique_ptr<MaterialProperty> prop) const{
    if (hasProperty(name))
        throw std::runtime_error("ERROR: A property of the same name already exists.");
    _props[name] = std::move(prop);
    _props[name]->setMaterial(this);
}

void Material::replaceProperty(Prop name, double val) noexcept{
    if (hasProperty(name)){
        _props[name] = std::make_unique<ConstantProperty>(val);
        _props[name]->setMaterial(this);
    }
}

void Material::replaceProperty(Prop name, std::function<double(const PropVars&)> func) noexcept{
    if (hasProperty(name)){
        _props[name] = std::make_unique<FunctionProperty>(func);
        _props[name]->setMaterial(this);
    }
}

void Material::replaceProperty(Prop name, std::unique_ptr<MaterialProperty> prop) noexcept{
    if (hasProperty(name)){
        _props[name].reset(prop.release());
        _props[name]->setMaterial(this);
    }
}

void Material::removeProperty(Prop name) noexcept{ _props.erase(name); }

double Material::computeProperty(Prop name, const PropVars& vars) const{
    if (!hasProperty(name))
        throw std::invalid_argument("ERROR: A property of this name has not been added.");
    return _props[name]->compute(vars);
}

double Material::rho(double P, double T) const noexcept{ return _eos->rho(P,T); }
double Material::drho_dP(double P, double T) const noexcept{ return _eos->drho_dP(P,T); }
double Material::drho_dT(double P, double T) const noexcept{ return _eos->drho_dT(P,T); }
double Material::beta(double P, double T) const noexcept{ return _eos->beta(P,T); }
double Material::Cp(double P, double T) const noexcept{ return _eos->Cp(P,T); }
double Material::dCp_dT(double P, double T) const noexcept{ return _eos->dCp_dT(P,T); }
double Material::k(double P, double T) const noexcept{ return _eos->k(P,T); }
double Material::mu(double P, double T) const noexcept{ return _eos->mu(P,T); }
double Material::Pr(double P, double T) const noexcept{ return _eos->Pr(P,T); }
double Material::H(double P, double T) const noexcept{ return _eos->H(P,T); }