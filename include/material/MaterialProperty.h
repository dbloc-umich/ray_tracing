#ifndef MATERIAL_PROPERTY_H
#define MATERIAL_PROPERTY_H

#include <vector>

class MaterialProperty{
    public:
    MaterialProperty() {};
    MaterialProperty(const MaterialProperty&) = delete;
    MaterialProperty(MaterialProperty&&) = default;
    virtual ~MaterialProperty() = default;
    MaterialProperty& operator=(const MaterialProperty&) = delete;
    MaterialProperty& operator=(MaterialProperty&&) = default;

    virtual double compute(const std::vector<double>& = {}) const = 0;
};

#endif