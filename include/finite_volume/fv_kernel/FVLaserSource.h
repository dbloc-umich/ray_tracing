#ifndef FV_LASER_SOURCE_H
#define FV_LASER_SOURCE_H

#include "FVKernel.h"

class Ray;
class Shape;
class FVLaserSource: public FVKernel{
    public:
    FVLaserSource(std::shared_ptr<Material> mat, Eigen::Index s, Ray& ray, Shape* shape);
    Eigen::VectorXd computeResidual(const FVStateMesh& u) const override;

    protected:
    Ray& _ray;
    Shape* _shape;
    double _I0;
};

#endif