#ifndef MATERIAL_H
#define MATERIAL_H

#include <map>
#include <string>
#include <vector>

class complex;
class Material{
    public:
    Material() {};

    protected:
    static const std::vector<std::string> validProperties;
    std::map<std::string, double> properties;
};

const auto Material::validProperties = {"extinction_coefficient",
                                        "refractive_index"};
#endif