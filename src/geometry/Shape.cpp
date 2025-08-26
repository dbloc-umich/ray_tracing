#include "Shape.h"
#include <exception>

Shape::Shape(double Sigma_t, double refrac):
    _Sigma_t(Sigma_t),
    _refrac(refrac)
    //_parent(nullptr)
{
    if (Sigma_t < 0.0) throw std::invalid_argument("ERROR: Cross section must be non-negative.");
    if (refrac < 1.0) throw std::invalid_argument("ERROR: Refractive index must be at least unity.");
}

void Shape::setSigma_t(double Sigma_t){
    if (Sigma_t < 0.0) throw std::invalid_argument("ERROR: Cross section must be non-negative.");
    _Sigma_t = Sigma_t;
}

void Shape::setRefractive(double refrac){
    if (refrac < 1.0) throw std::invalid_argument("ERROR: Refractive index must be at least unity.");
    _refrac = refrac;
}