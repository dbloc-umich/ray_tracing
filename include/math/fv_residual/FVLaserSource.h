#ifndef FV_LASER_SOURCE_H
#define FV_LASER_SOURCE_H

#include "FVResidual.h"

class Ray;
class Shape;
class FVLaserSource: public FVResidual{
    public:
    FVLaserSource(std::shared_ptr<Material> mat, Ray& ray, Shape* shape);
    Eigen::VectorXd computeResidual(const FVStateMesh& u) const override;

    protected:
    Ray& _ray;
    Shape* _shape;
    double _I0;
};

#endif