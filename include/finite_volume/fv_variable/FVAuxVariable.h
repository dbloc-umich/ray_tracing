#ifndef FV_AUX_VARIABLE_H
#define FV_AUX_VARIABLE_H

#include "FVVariable.h"
#include <memory>
class FVAuxKernel;
class FVAuxVariable: public FVVariable{
    public:
    FVAuxVariable(const std::string& name, unsigned short dim = 1): FVVariable(name, dim) {}
    bool aux() const noexcept override{ return true; }
    void setAuxKernel(std::unique_ptr<FVAuxKernel> auxKernel){ _auxKernel = std::move(auxKernel); }
    double computeValue() const;

    protected:
    std::unique_ptr<FVAuxKernel> _auxKernel;
};

#endif