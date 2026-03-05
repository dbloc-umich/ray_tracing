// #ifndef NONLINEAR_VARIABLE_H
// #define NONLINEAR_VARIABLE_H

// #include "Variable.h"
// #include "Eigen/Dense"
// #include <memory>
// #include <vector>

// class Kernel;
// class StateMesh;
// class NonlinearVariable: public Variable{
//     public:
//     NonlinearVariable(const std::string& name, unsigned short dim = 1): Variable(name, dim) {}
//     bool aux() const noexcept override{ return false; }
//     void addKernel(std::unique_ptr<Kernel>);
//     Eigen::VectorXd computeResidual(const StateMesh& u) const;

//     protected:
//     std::vector<std::unique_ptr<Kernel>> _kernels;
// };

// #endif