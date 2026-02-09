#ifndef CONSTANT_PROPERTY_H
#define CONSTANT_PROPERTY_H

#include "MaterialProperty.h"

class ConstantProperty : public MaterialProperty{
    public:
    ConstantProperty(double val): _val(val) {}
    ConstantProperty(const ConstantProperty&) = delete;
    ConstantProperty(ConstantProperty&& other): _val(other._val) {}
    ConstantProperty& operator=(const ConstantProperty&) = delete;
    ConstantProperty& operator=(ConstantProperty&& other){ _val = other._val; return *this; }

    double compute(const std::map<PropVariable, double>& = {}) const override{ return _val; };

    protected:
    double _val;
};

#endif