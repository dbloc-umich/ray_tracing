// An auxiliary kernel computes a quantity of interest given a state vector within a cell
#ifndef AUXKERNEL_H
#define AUXKERNEL_H

#include "StateMesh.h"
class AuxKernel{
    public:
    AuxKernel() {}
    virtual ~AuxKernel() = default;
    // Best if the map is returned by StateMesh::matProp
    virtual double computeValue(const std::map<std::string, double>& u) const = 0;
};

#endif