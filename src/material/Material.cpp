#include "Material.h"
#include "MaterialProperty.h"
#include <algorithm>
#include <stdexcept>

const std::set<std::string> Material::_validProps{"extinction_coefficient",
                                                  "refractive_index"};

Material::Material() = default;
Material::~Material() = default;
Material::Material(Material&&) = default;
Material& Material::operator=(Material&&) = default;

void Material::addProperty(std::string name, std::unique_ptr<MaterialProperty> prop) const{
    if (std::find(_validProps.cbegin(), _validProps.cend(), name) == _validProps.cend())
        throw std::invalid_argument("ERROR: Invalid property name.");

    if (_props.find(name) != _props.cend())
        throw std::runtime_error("ERROR: A property of the same name already exists. Consider using replaceProperty.");
    _props[name] = std::move(prop);
}

void Material::replaceProperty(std::string name, std::unique_ptr<MaterialProperty>& prop){
    if (_props.find(name) != _props.cend()) _props[name].swap(prop);
    prop.release();
}

void Material::removeProperty(std::string name){
    _props.erase(name);
}

double Material::computeProperty(std::string name, const std::vector<double>& vars) const{
    if (_props.find(name) == _props.cend())
        throw std::invalid_argument("ERROR: A property of this name has not been added.");
    return _props[name]->compute(vars);
}