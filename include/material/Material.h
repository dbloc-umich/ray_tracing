#ifndef MATERIAL_H
#define MATERIAL_H

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

class MaterialProperty;
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
    
    void addProperty(std::string name, std::unique_ptr<MaterialProperty> prop) const;
    void replaceProperty(std::string name, std::unique_ptr<MaterialProperty>& prop);
    void removeProperty(std::string name);
    double computeProperty(std::string name, const std::vector<double>& vars = {}) const;

    protected:
    mutable std::map<std::string, std::unique_ptr<MaterialProperty>> _props;
    static const std::set<std::string> _validProps;
};

#endif