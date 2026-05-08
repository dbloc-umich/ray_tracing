#ifndef EXPLICIT_TIME_INTEGRATOR_H
#define EXPLICIT_TIME_INTEGRATOR_H

#include "TimeIntegrator.h"
class ExplicitTimeIntegrator: public TimeIntegrator{
    public:
    virtual ~ExplicitTimeIntegrator() = default;  
};

#endif