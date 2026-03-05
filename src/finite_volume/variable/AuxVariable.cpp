// #include "AuxVariable.h"
// #include "AuxKernel.h"

// void AuxVariable::setAuxKernel(std::unique_ptr<AuxKernel> auxKernel){
//     if (auxKernel) _auxKernel = std::move(auxKernel);
// }

// double AuxVariable::computeValue() const{
//     if (_auxKernel) return _auxKernel->computeValue();
//     return 0;
// }