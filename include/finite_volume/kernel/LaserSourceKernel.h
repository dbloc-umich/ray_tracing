#ifndef LASER_SOURCE_KERNEL_H
#define LASER_SOURCE_KERNEL_H

#include "Kernel.h"

class Ray;
class Shape;
class LaserSourceKernel: public Kernel{
    public:
    LaserSourceKernel(std::shared_ptr<Material> mat, std::shared_ptr<Material> eMat, Ray& ray, Shape* shape);
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;
    Eigen::VectorXi stateID(const StateMesh& u) const noexcept override;

    protected:
    std::shared_ptr<Material> _eMat;
    Ray& _ray;
    Shape* _shape;
    double _I0; // initial laser intensity
};

#endif