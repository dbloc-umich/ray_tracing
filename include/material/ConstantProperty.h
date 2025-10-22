#ifndef CONSTANT_PROPERTY_H
#define CONSTANT_PROPERTY_H

#include "MaterialProperty.h"

class ConstantProperty : public MaterialProperty{
    public:
    ConstantProperty(double val): _val(val) {}
    double compute(const std::vector<double>& = {}) const override{ return _val; };

    protected:
    double _val;
};

#endif