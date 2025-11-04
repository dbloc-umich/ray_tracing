#include "Ray.h"
#include "Constants.h"
#include "Shape.h"

#define _USE_MATH_DEFINES
#include <cmath>

Ray::Ray(Eigen::Vector3d p, UnitVector3d dir, double I, double lambda, Shape* host):
    _p(p),
    _dir(dir),
    _I(I),
    _lambda(lambda),
    _host(host)
{
    if (I < 0.0) throw std::invalid_argument("ERROR: Negative intensity.");
    if (lambda <= 0.0) throw std::invalid_argument("ERROR: Non-positive wavelength.");
    if (host && !(host->encloses(p) || host->surfaceContains(p)))
        throw std::invalid_argument("ERROR: Shape pointer does not point to a Shape that contains the starting point.");
}

void Ray::setIntensity(double I){
    if (I < 0.0) throw std::invalid_argument("ERROR: Negative intensity.");
    _I = I;
}

double Ray::frequency() const noexcept{
    // angular frequency in units of radians/s
    return 2.0 * M_PI * constants::c / _lambda;
}

void Ray::setWavelength(double lambda){
    if (lambda <= 0.0) throw std::invalid_argument("ERROR: Non-positive wavelength.");
    _lambda = lambda;
}

void Ray::setHost(Shape* host){
    if (host && !(host->encloses(_p) || host->surfaceContains(_p)))
        throw std::invalid_argument("ERROR: Shape pointer does not point to a Shape that contains the starting point.");
    _host = host; 
}