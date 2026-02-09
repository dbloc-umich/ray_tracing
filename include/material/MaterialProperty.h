#ifndef MATERIAL_PROPERTY_H
#define MATERIAL_PROPERTY_H

#include <map>

class Material;
enum class PropVariable;
class MaterialProperty{
    public:
    MaterialProperty() {};
    MaterialProperty(const MaterialProperty&) = delete;
    MaterialProperty(MaterialProperty&&){};
    virtual ~MaterialProperty() = default;
    MaterialProperty& operator=(const MaterialProperty&) = delete;
    MaterialProperty& operator=(MaterialProperty&&){ return *this; }

    void setMaterial(const Material* mat) noexcept{ _mat = mat; }
    virtual double compute(const std::map<PropVariable, double>& = {}) const = 0;

    protected:
    const Material* _mat; // pointer to a const "parent" Material object
};

#endif