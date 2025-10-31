#ifndef FUNCTION_PROPERTY_H
#define FUNCTION_PROPERTY_H

#include "MaterialProperty.h"
#include <functional>

class FunctionProperty : public MaterialProperty{
    public:
    FunctionProperty(std::function<double(const std::vector<double>&)> f): _f(f) {};
    double compute(const std::vector<double>& var = {}) const override{ return _f(var); };

    protected:
    std::function<double(const std::vector<double>&)> _f;
};

#endif