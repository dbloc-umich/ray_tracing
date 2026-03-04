#ifndef FV_AUX_KERNEL_H
#define FV_AUX_KERNEL_H

class FVParameters;
class FVAuxKernel{
    public:
    FVAuxKernel(const FVParameters& param) {}
    virtual double computeValue() = 0;
};

#endif