#include "RayTracing.h"
#include "Constants.h"
#include "Material.h"
#include "Ray.h"
#include "Shape.h"
#include "Tree.h"

#include <random>
#include <stack>

//#define MONITOR

namespace{
    static constexpr double IThreshold = 1e-9;
    static std::default_random_engine rng(64);
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    template<typename T>
    Shape* attenuate(const Tree<T>& tree, Ray& ray, double& s){
        // A single attenuation event, update the position and intensity of the current ray
        // Returns the next host and update the distance to that host

        Shape* next = tree.nextShape(ray.position(), ray.direction(), ray.host(), s);
#ifdef MONITOR
    std::cout << "The ray intensity is " << ray.intensity() << ", ";
    if (!ray.host()){
        if (!next) std::cout << "Point " << ray.position() << " traveling at " << ray.direction() << " does not reach any other Shape.";
        else std::cout << "Point " << ray.position() << " traveling at " << ray.direction() << " reaches " << *next << " after a distance of " << s << ".";
    } else{
        if (!next) std::cout << "Point " << ray.position() << " on " << *ray.host() << " traveling at " << ray.direction() << " does not reach any other Shape.";
        else std::cout << "Point " << ray.position() << " on " << *ray.host() << " traveling at " << ray.direction() << " reaches " << *next << " after a distance of " << s << ".";
    }
    std::cout << std::endl;
#endif
        if (next){
            ray.setPoistion(ray.position().advance(ray.direction(), s));
            if (ray.host() == next){
                double alpha = 0.0;
                if (next->hasProperty(Prop::attenuationCoefficient)){
                    alpha = next->computeProperty(Prop::attenuationCoefficient, {ray.wavelength()});
                } else if (next->hasProperty(Prop::extinctionCoefficient)){
                    alpha = next->computeProperty(Prop::extinctionCoefficient, {ray.wavelength()});
                    alpha *= 2.0 * ray.frequency() / constants::c;
                }
                ray.setIntensity(ray.intensity()*exp(-alpha*s)); // moving within a Shape, implicit absorption;
            }
        }
        return next;
    }

    template<typename T>
    std::vector<Ray> penetrate(const Tree<T>& tree, Ray& ray, Shape* next, double s, const HopMode& mode){
        /**
         * A single penetration event at the interface between two medium
         * If an empty vector is returned, the ray only takes the one path updated in ray
         * If a vector of one Ray is returned, the ray is split into two paths,
         *  one stored in ray and one returned in the vector
        **/

        double n1, n2;
        if (!ray.host()){ // current points to empty space
            n1 = 1.0;
            n2 = next->computeProperty(Prop::refractiveIndex, {ray.wavelength()});
        } else if (s == 0.0){
            // particle enters directly from current to next
            n1 = ray.host()->computeProperty(Prop::refractiveIndex, {ray.wavelength()});
            n2 = next->computeProperty(Prop::refractiveIndex, {ray.wavelength()});
        } else if (ray.host() != next){
            // particle travels from one Node to another with vacuum in between them
            n1 = 1.0;
            n2 = next->computeProperty(Prop::refractiveIndex, {ray.wavelength()});
        } else{
            // particle travels within the same Node
            n1 = ray.host()->computeProperty(Prop::refractiveIndex, {ray.wavelength()});
            auto pseudoNext = tree.nextShape(ray.position(), ray.direction(), ray.host(), s);
            n2 = (s == 0) ? pseudoNext->computeProperty(Prop::refractiveIndex, {ray.wavelength()}) : 1.0;
        }

        Direction normal = next->normal(ray.position());
        Direction transmit = refracted(ray.direction(), normal, n1, n2);
        Direction reflect = reflected(ray.direction(), normal);
        
        ray.setDirection(reflect);
        if (!transmit) return {}; // Total internal reflection    
        
        // Transmission and reflection treatment
        double cosi = fabs(reflect.dot(normal));
        double cost = fabs(transmit.dot(normal));
        double Rs = (n1*cosi - n2*cost)/(n1*cosi + n2*cost);
        double Rp = (n1*cost - n2*cosi)/(n1*cost + n2*cosi);
        double R = 0.5*(Rs*Rs + Rp*Rp);
        if (mode == HopMode::Roulette){
            if (dist(rng) < R) ray.setDirection(reflect);
            else ray.setDirection(transmit);
            return {};
        }
        std::vector<Ray> ray2{ray};
        ray.setDirection(transmit);
        ray.setIntensity(ray.intensity()*(1-R));
        ray2[0].setDirection(reflect);
        ray2[0].setIntensity(ray2[0].intensity()*R);
        return ray2;
    }
}

template<typename T>
double intensity(const Tree<T>& tree, Ray& ray, const HopMode& mode) noexcept{
    if (!tree) return ray.intensity();
    if (ray.intensity() == 0.0) return 0.0;
    
    double I0 = ray.intensity();
    if (mode != HopMode::Detailed){
        double s;
        Shape* next = attenuate(tree, ray, s);
        while(next){
            if (ray.intensity() < IThreshold*I0){
                if (dist(rng) < ray.intensity()/I0) ray.setIntensity(IThreshold*I0);
                else return 0.0;
            }
            if (mode == HopMode::Roulette) penetrate(tree, ray, next, s, mode);
            ray.setHost(next);
            next = attenuate(tree, ray, s);
        }
        return ray.intensity();
    }

    std::stack<Ray> stack;
    stack.push(std::move(ray));
    double sum = 0.0;
    while (!stack.empty()){
        ray = stack.top();
        stack.pop();

        if (ray.intensity() < IThreshold*I0){
            if (dist(rng) <= ray.intensity()/I0) ray.setIntensity(IThreshold*I0);
            else continue;
        }

        double s; // placeholder
        Shape* next = attenuate(tree, ray, s);
        if (!next) sum += ray.intensity();
        else{
            auto ray2 = penetrate(tree, ray, next, s, mode);
            ray.setHost(next);
            stack.push(std::move(ray));
            if (!ray2.empty()){
                ray2[0].setHost(next);
                stack.push(std::move(ray2[0]));
            }
        }
    }
    return sum;
}

// Explicit instantiations
#include "Node.h"
template double intensity(const Tree<Node>&, Ray&, const HopMode&);

#include "UNode.h"
template double intensity(const Tree<UNode>&, Ray&, const HopMode&);