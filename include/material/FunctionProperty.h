#ifndef FUNCTION_PROPERTY_H
#define FUNCTION_PROPERTY_H

#include "MaterialProperty.h"
#include <functional>

class FunctionProperty : public MaterialProperty{
    public:
    FunctionProperty(std::function<double(const std::map<PropVariable, double>&)> f): _f(f) {};
    FunctionProperty(const FunctionProperty&) = delete;
    FunctionProperty(FunctionProperty&& other): _f(std::move(other._f)) {}
    FunctionProperty& operator=(const FunctionProperty&) = delete;
    FunctionProperty& operator=(FunctionProperty&& other){ _f = std::move(other._f); return *this; }

    double compute(const std::map<PropVariable, double>& var = {}) const override{ return _f(var); };

    protected:
    std::function<double(const std::map<PropVariable, double>&)> _f;
};

#endif