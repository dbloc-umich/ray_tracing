#ifndef LASER_SOURCE_KERNEL_H
#define LASER_SOURCE_KERNEL_H

#include "Kernel.h"

class Ray;
class Shape;
class LaserSourceKernel: public Kernel{
    public:
    LaserSourceKernel(std::shared_ptr<Material> mat, Ray& ray, Shape* shape);
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;

    protected:
    Ray& _ray;
    Shape* _shape;
    double _I0; // initial laser intensity
};

#endif