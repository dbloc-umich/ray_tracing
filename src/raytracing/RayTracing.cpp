#include "RayTracing.h"
#include "Direction.h"
#include "Point.h"
#include "Shape.h"
#include "Tree.h"

#include <cmath>
#include <random>

namespace{
    static constexpr double IThreshold = 1e-9;
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
        return ortho*sine + sgn*normal*std::sqrt(1-sine*sine);
    }
    return Direction(nullptr);
}

template<typename T>
double intensity(const Tree<T>& tree, Point p, const Direction& dir, bool isRefracted, bool isReflected, Shape* current, double initial){
    if (initial == 0.0) return 0.0;
    if (initial < IThreshold){
        if (dist(rng) <= initial/IThreshold) return intensity(tree, p, dir, isRefracted, isReflected, current, IThreshold);
        return 0.0;
    }

    double s; // placeholder
    Shape* next = tree.nextShape(p, dir, current, s);
    // if (!current){
    //     if (!next) std::cout << "Point " << p << " traveling at " << dir << " does not reach any other Shape." << std::endl;
    //     else std::cout << "Point " << p << " traveling at " << dir << " reaches " << *next << " after a distance of " << s << "." << std::endl;
    // } else{
    //     if (!next) std::cout << "Point " << p << " on " << *current << " traveling at " << dir << " does not reach any other Shape." << std::endl;
    //     else std::cout << "Point " << p << " on " << *current << " traveling at " << dir << " reaches " << *next << " after a distance of " << s << "." << std::endl;
    // }
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
        auto pseudoNext = tree.nextShape(p, dir, current, s);
        n2 = (s == 0) ? pseudoNext->refractive() : 1.0;
    }

    Direction normal = next->normal(p);
    Direction refract = refracted(dir, normal, n1, n2);
    Direction reflect = reflected(dir, normal);

    if (isRefracted && refract){
        // if the refration is calculated and the refracted ray exists
        if (isReflected){
            double cosi = fabs(dir.dot(normal));
            double cost = fabs(refract.dot(normal));
            double Rs = (n1*cosi - n2*cost)/(n1*cosi + n2*cost);
            double Rp = (n1*cost - n2*cosi)/(n1*cost + n2*cosi);
            double R = 0.5*(Rs*Rs + Rp*Rp);
            return intensity(tree, p, refract, isRefracted, isReflected, next, initial*(1-R))
                 + intensity(tree, p, reflect, isRefracted, isReflected, next, initial*R);
        }
        return intensity(tree, p, refract, isRefracted, isReflected, next, initial);
    }
    // only reflection, no refraction
    return intensity(tree, p, reflect, isRefracted, isReflected, next, initial);
}

// Explicit instantiations
#include "Node.h"
template double intensity<Node>(const Tree<Node>&, Point, const Direction&, bool, bool, Shape*, double);

#include "UNode.h"
template double intensity<UNode>(const Tree<UNode>&, Point, const Direction&, bool, bool, Shape*, double);
