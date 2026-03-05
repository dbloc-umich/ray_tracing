// #ifndef AUX_VARIABLE_H
// #define AUX_VARIABLE_H

// #include "Variable.h"
// #include <memory>
// class AuxKernel;
// class AuxVariable: public Variable{
//     public:
//     AuxVariable(const std::string& name, unsigned short dim = 1): Variable(name, dim) {}
//     bool aux() const noexcept override{ return true; }
//     void setAuxKernel(std::unique_ptr<AuxKernel> auxKernel);
//     double computeValue() const;

//     protected:
//     std::unique_ptr<AuxKernel> _auxKernel;
// };

// #endif