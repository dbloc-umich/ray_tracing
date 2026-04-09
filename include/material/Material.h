#ifndef MATERIAL_H
#define MATERIAL_H

#include <functional>
#include <map>
#include <memory>

class MaterialProperty; // for all other properties
class Material{
    protected:
    using PropVars = std::map<std::string, double>;

    public:
    Material();
    Material(const Material&) = delete;
    Material(Material&&);
    ~Material();
    Material& operator=(const Material&) = delete;
    Material& operator=(Material&&);
    // the unique_ptr does not have access to a default deleter of an incomplete class yet,
    // so any automatically generated functions cannot be defined until the class is complete.

    bool hasProperty(const std::string& name) const noexcept;

    void addProperty(const std::string& name) const;
    void addProperty(const std::string& name, double val) const;
    void addProperty(const std::string& name, std::function<double(const PropVars&)> func) const;
    void addProperty(const std::string& name, std::unique_ptr<MaterialProperty> prop) const;

    void replaceProperty(const std::string& name, double val) noexcept;
    void replaceProperty(const std::string& name, std::function<double(const PropVars&)> func) noexcept;
    void replaceProperty(const std::string& name, std::unique_ptr<MaterialProperty> prop) noexcept;

    void removeProperty(const std::string& name) noexcept;
    virtual double computeProperty(const std::string& name, const PropVars& vars = {}) const;

    // // Functions from eos
    // virtual double T_from_H(double H) const = 0;

    protected:
    mutable std::map<std::string, std::unique_ptr<MaterialProperty>> _props;
};

#endif