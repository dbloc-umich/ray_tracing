#include "Material.h"
#include "MaterialProperty.h"
#include <algorithm>
#include <stdexcept>

Material::Material() = default;
Material::~Material() = default;
Material::Material(Material&&) = default;
Material& Material::operator=(Material&&) = default;

void Material::addProperty(Prop name, std::unique_ptr<MaterialProperty> prop) const{
    if (_props.find(name) != _props.cend())
        throw std::runtime_error("ERROR: A property of the same name already exists. Consider using replaceProperty.");
    _props[name] = std::move(prop);
}

void Material::replaceProperty(Prop name, std::unique_ptr<MaterialProperty>& prop) noexcept{
    if (_props.find(name) != _props.cend()) _props[name].swap(prop);
    prop.release();
}

void Material::removeProperty(Prop name) noexcept{
    _props.erase(name);
}

double Material::computeProperty(Prop name, const std::vector<double>& vars) const{
    if (_props.find(name) == _props.cend())
        throw std::invalid_argument("ERROR: A property of this name has not been added.");
    return _props[name]->compute(vars);
}