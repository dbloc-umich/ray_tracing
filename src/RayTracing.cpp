#include "RayTracing.h"
#include "Direction.h"
#include "kDTree.h"

#include <cmath>
#include <random>

namespace{
    static constexpr double IThreshold = 1e-6;
    static std::default_random_engine rng(64);
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
}

Direction reflected(const Direction& in, const Direction& normal){ return in - 2*(normal.dot(in))*normal; }
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2){
    if (in.isParallel(normal) || in.isOrthogonal(normal)) return in;
    if (n1 == n2) return in;

    Vector c = in.cross(normal)*(n1/n2);
    double sine = c.norm();
    if (sine <= 1.0){
        Direction ortho = normal.cross(c); // component of the refracted vector that is orthogonal to normal
        int sgn = (in.dot(normal) > 0.0) ? 1 : -1;
        return ortho*sine + sgn*normal*sqrt(1-sine*sine);
    }
    return Direction(nullptr);
}

double intensity(const kDTree& tree, Point p, const Direction& dir, Shape* current, double initial){
    if (initial < IThreshold){
        double prob = initial/IThreshold; // survival probability
        std::cout << prob << std::endl;
        if (dist(rng) <= prob) return intensity(tree, p, dir, nullptr, IThreshold);
        return 0.0;
    }
    double s; // placeholder
    Shape* next = tree.nextNode(p, dir, current, s);
    if (!next) return initial;
    {
        //std::cout << "Recovering " << initial << std::endl;
    }

    p.advance(dir, s);
    if (current == next) initial *= exp(-next->Sigma_t()*s); // moving within a Shape, implicit absorption

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
    Direction refract = refracted(dir, normal, n1, n2);
    Direction reflect = reflected(dir, normal);
    if (refract){
        double cosi = fabs(dir.dot(normal));
        double cost = fabs(refract.dot(normal));
        double R = (n1*cosi - n2*cost)/(n1*cosi + n2*cost);
        R *= R; // reflection coefficient
        return intensity(tree, p, reflect, next, initial*R) + intensity(tree, p, refract, next, initial*(1-R));
    }
    return intensity(tree, p, reflect, next, initial); // total internal reflection
}