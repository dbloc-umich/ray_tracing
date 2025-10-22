#include "Ray.h"
#include "Shape.h"

Ray::Ray(Point p, Direction dir, double I, double lambda, Shape* host):
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

void Ray::setWavelength(double lambda){
    if (lambda <= 0.0) throw std::invalid_argument("ERROR: Non-positive wavelength.");
    _lambda = lambda;
}

void Ray::setHost(Shape* host){
    if (host && !(host->encloses(_p) || host->surfaceContains(_p)))
        throw std::invalid_argument("ERROR: Shape pointer does not point to a Shape that contains the starting point.");
    _host = host; 
}