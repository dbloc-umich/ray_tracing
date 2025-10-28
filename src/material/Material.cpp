#include "Material.h"
#include "MaterialProperty.h"
#include <algorithm>
#include <stdexcept>

Material::Material() = default;
Material::~Material() = default;
Material::Material(Material&&) = default;
Material& Material::operator=(Material&&) = default;

bool Material::hasProperty(Prop name) const noexcept{
    return _props.count(name);
}

void Material::addProperty(Prop name, std::unique_ptr<MaterialProperty> prop) const{
    if (hasProperty(name))
        throw std::runtime_error("ERROR: A property of the same name already exists.");
    _props[name] = std::move(prop);
}

void Material::replaceProperty(Prop name, std::unique_ptr<MaterialProperty> prop) noexcept{
    if (hasProperty(name)) _props[name].reset(prop.release());
    else _props[name].swap(prop);
}

void Material::removeProperty(Prop name) noexcept{
    _props.erase(name);
}

double Material::computeProperty(Prop name, const std::vector<double>& vars) const{
    if (!hasProperty(name))
        throw std::invalid_argument("ERROR: A property of this name has not been added.");
    return _props[name]->compute(vars);
}