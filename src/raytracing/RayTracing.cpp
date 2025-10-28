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
}

template<typename T>
double intensity(const Tree<T>& tree, const Ray& ray, bool isRefracted, bool isReflected) noexcept{
    if (!tree) return ray.intensity();
    if (ray.intensity() == 0.0) return 0.0;
    
    std::stack<Ray> stack;
    stack.push(ray);
    double sum = 0.0;

    while (!stack.empty()){
        Ray top = stack.top();
        stack.pop();
        if (top.intensity() < IThreshold){
            if (dist(rng) <= top.intensity()/ray.intensity()) top.setIntensity(IThreshold*ray.intensity());
            else continue;
        }

        double s; // placeholder
        Shape* next = tree.nextShape(top.position(), top.direction(), top.host(), s);
#ifdef MONITOR
    std::cout << "The ray intensity is " << top.intensity() << ", ";
    if (!ray.host()){
        if (!next) std::cout << "Point " << top.position() << " traveling at " << top.direction() << " does not reach any other Shape.";
        else std::cout << "Point " << top.position() << " traveling at " << top.direction() << " reaches " << *next << " after a distance of " << s << ".";
    } else{
        if (!next) std::cout << "Point " << top.position() << " on " << *ray.host() << " traveling at " << top.direction() << " does not reach any other Shape.";
        else std::cout << "Point " << top.position() << " on " << *ray.host() << " traveling at " << top.direction() << " reaches " << *next << " after a distance of " << s << ".";
    }
    std::cout << std::endl;
#endif
        if (!next) sum += top.intensity();
        else{
            top.setPoistion(top.position().advance(top.direction(), s));
            if (top.host() == next){
                double alpha = 0.0;
                if (next->hasProperty(Prop::attenuationCoefficient)){
                    alpha = next->computeProperty(Prop::attenuationCoefficient, {top.wavelength()});
                } else if (next->hasProperty(Prop::extinctionCoefficient)){
                    alpha = next->computeProperty(Prop::extinctionCoefficient, {top.wavelength()});
                    alpha *= 2.0 * top.frequency() / constants::c;
                }
                top.setIntensity(top.intensity()*exp(-alpha*s)); // moving within a Shape, implicit absorption;
            }
            if (!isRefracted && !isReflected){
                top.setHost(next);
                stack.push(std::move(top));
            } else{
                double n1, n2;
                if (!top.host()){ // current points to empty space
                    n1 = 1.0;
                    n2 = next->computeProperty(Prop::refractiveIndex, {top.wavelength()});
                } else if (s == 0.0){
                    // particle enters directly from current to next
                    n1 = top.host()->computeProperty(Prop::refractiveIndex, {top.wavelength()});
                    n2 = next->computeProperty(Prop::refractiveIndex, {top.wavelength()});
                } else if (top.host() != next){
                    // particle travels from one Node to another with vacuum in between them
                    n1 = 1.0;
                    n2 = next->computeProperty(Prop::refractiveIndex, {top.wavelength()});
                } else{
                    // particle travels within the same Node
                    n1 = top.host()->computeProperty(Prop::refractiveIndex, {top.wavelength()});
                    auto pseudoNext = tree.nextShape(top.position(), top.direction(), top.host(), s);
                    n2 = (s == 0) ? pseudoNext->computeProperty(Prop::refractiveIndex, {top.wavelength()}) : 1.0;
                }

                Direction normal = next->normal(top.position());
                Direction refract = refracted(top.direction(), normal, n1, n2);
                Direction reflect = reflected(top.direction(), normal);

                top.setHost(next);
                if (isRefracted && refract){
                    // if the refration is calculated and the refracted ray exists
                    if (isReflected){
                        double cosi = fabs(top.direction().dot(normal));
                        double cost = fabs(refract.dot(normal));
                        double Rs = (n1*cosi - n2*cost)/(n1*cosi + n2*cost);
                        double Rp = (n1*cost - n2*cosi)/(n1*cost + n2*cosi);
                        double R = 0.5*(Rs*Rs + Rp*Rp);

                        Ray reflectedRay = top;
                        reflectedRay.setDirection(reflect);
                        reflectedRay.setIntensity(top.intensity()*R);
                        stack.push(std::move(reflectedRay));

                        top.setDirection(refract); // top is now the reflected beam;
                        top.setIntensity(top.intensity()*(1-R));
                        stack.push(std::move(top)); // refracted ray
                    } else{
                        // only refraction, no reflection
                        top.setDirection(refract);
                        stack.push(std::move(top));
                    }
                } else{
                    // only reflection, no refraction
                    top.setDirection(reflect);
                    stack.push(std::move(top));
                }
            }
        }
    }
    return sum;
}

// Explicit instantiations
#include "Node.h"
template double intensity(const Tree<Node>&, const Ray&, bool, bool);

#include "UNode.h"
template double intensity(const Tree<UNode>&, const Ray&, bool, bool);