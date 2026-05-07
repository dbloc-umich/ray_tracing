#ifndef ELECTRON_MATERIAL_H
#define ELECTRON_MATERIAL_H

#include "Material.h"
#include "ElectronEquationOfState.h"
#include <string>

class ElectronMaterial: public Material{
    public:
    ElectronMaterial(const EquationOfState&);
    double computeProperty(const std::string& name, const PropVars& vars = {}) const override;
    const EquationOfState& eos() const noexcept override;

    protected:
    static const ElectronEquationOfState _eos;
    const EquationOfState& _nEOS; // eos of the neutral species
};
#endif