#ifndef LASER_SOURCE_H
#define LASER_SOURCE_H

#include "Kernel.h"

class Ray;
class Shape;
class LaserSource: public Kernel{
    public:
    LaserSource(std::shared_ptr<Material> mat, const Eigen::VectorXi& s, Ray& ray, Shape* shape);
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;

    protected:
    Ray& _ray;
    Shape* _shape;
    double _I0; // initial laser intensity
};

#endif