#include "RayTracing.h"
#include "Direction.h"
#include "Octree.h"

#include <cmath>
#include <random>

namespace{
    static constexpr double IThreshold = 1e-6;
    static std::default_random_engine rng(64);
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
}

Direction reflected(const Direction& in, const Direction& normal){ return in - 2*(normal.dot(in))*normal; }
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2){
    if (in.isParallel(normal) || in.isOrthogonal(normal) || n1 == n2) return in;

    Vector c = in.cross(normal)*(n1/n2);
    double sine = c.norm();
    if (sine <= 1.0){
        Direction ortho = normal.cross(c); // component of the refracted vector that is orthogonal to normal
        int sgn = in.dot(normal) > 0.0 ? 1 : -1;
        return ortho*sine + sgn*normal*sqrt(1-sine*sine);
    }
    return Direction(nullptr);
}

double intensity(const Octree& tree, Point p, const Direction& dir, bool isRefracted, bool isReflected, Shape* current, double initial){
    if (initial < IThreshold){
        if (dist(rng) <= initial/IThreshold) return intensity(tree, p, dir, current, IThreshold);
        return 0.0;
    }
    double s; // placeholder
    Shape* next = tree.nextNode(p, dir, current, s);
    if (!next) return initial;

    p.advance(dir, s);
    if (current == next) initial *= exp(-next->Sigma_t()*s); // moving within a Shape, implicit absorption
    if (!isRefracted && !isReflected) return intensity(tree, p, dir, isRefracted, isReflected, next, initial);

    double n1, n2;
    if (!current){ // current points to empty space
        n1 = 1.0;
        n2 = next->refractive();
    } else if (s == 0.0){
        // particle enters directly from current to next
        n1 = current->refractive();
        n2 = next->refractive();
    } else if (current != next){
        // particle travels from one Node to another with vacuum in between them
        n1 = 1.0;
        n2 = next->refractive();
    } else{
        // particle travels within the same Node
        n1 = current->refractive();
        auto pseudoNext = tree.nextNode(p, dir, current, s);
        n2 = (s == 0) ? pseudoNext->refractive() : 1.0;
    }

    Direction normal = next->normal(p);
    Direction refract = isRefracted ? refracted(dir, normal, n1, n2) : Direction(nullptr);
    Direction reflect = isReflected ? reflected(dir, normal) : Direction(nullptr);

    if (isRefracted && refract){
        // if the refration is calculated and the refracted ray exists
        double cosi = fabs(dir.dot(normal));
        double cost = fabs(refract.dot(normal));
        double Rs = (n1*cosi - n2*cost)/(n1*cosi + n2*cost);
        double Rp = (n1*cost - n2*cosi)/(n1*cost + n2*cosi);
        double R = 0.5*(Rs*Rs + Rp*Rp); // reflectance
        double refractedIntensity = intensity(tree, p, refract, isRefracted, isReflected, next, initial*(1-R));
        double reflectedIntensity = isReflected ? intensity(tree, p, reflect, isRefracted, isReflected, next, initial*R) : 0.0;
        return refractedIntensity + reflectedIntensity;
    }
    // only reflection, no refraction
    return intensity(tree, p, reflect, isRefracted, isReflected, next, initial);
}