#ifndef INTERPOLATED_PROPERTY_H
#define INTERPOLATED_PROPERTY_H

#include "MaterialProperty.h"

#include <string>

class InterpolatedProperty : public MaterialProperty{
    public:
    InterpolatedProperty(std::string file);
    double compute(const std::vector<double>& = {}) const override;

    protected:
    std::vector<double> _grid;
    std::vector<double> _val;

    void loadFile(std::string file);
};

#endif