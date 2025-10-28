#ifndef MATERIAL_H
#define MATERIAL_H

#include <map>
#include <memory>
#include <vector>

class MaterialProperty;
enum class Prop{attenuationCoefficient,
                extinctionCoefficient,
                refractiveIndex};

class Material{
    public:
    Material();
    Material(const Material&) = delete;
    Material(Material&&);
    ~Material();
    Material& operator=(const Material&) = delete;
    Material& operator=(Material&&);
    // the unique_ptr does not have access to a default deleter of an incomplete class yet,
    // so any automatically generated functions cannot be defined until the class is complete.
    
    bool hasProperty(Prop name) const noexcept;
    void addProperty(Prop name, std::unique_ptr<MaterialProperty> prop) const;
    void replaceProperty(Prop name, std::unique_ptr<MaterialProperty> prop) noexcept;
    void removeProperty(Prop name) noexcept;
    double computeProperty(Prop name, const std::vector<double>& vars = {}) const;

    protected:
    mutable std::map<Prop, std::unique_ptr<MaterialProperty>> _props;
};

#endif