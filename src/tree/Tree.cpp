#include "Tree.h"

#include <cmath>
#include <limits>
#include <random>
#include <stack>

#include "Box.h"
#include "Point.h"
#include "Ray.h"

//#define MONITOR

namespace{
    static constexpr double IThreshold = 1e-9;
    static std::default_random_engine rng(64);
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
}

template <class T>
double Tree<T>::xMin() const noexcept{ return _root ? _root->xMin() : NAN; }

template <class T>
double Tree<T>::xMax() const noexcept{ return _root ? _root->xMax() : NAN; }

template <class T>
double Tree<T>::yMin() const noexcept{ return _root ? _root->yMin() : NAN; }

template <class T>
double Tree<T>::yMax() const noexcept{ return _root ? _root->yMax() : NAN; }

template <class T>
double Tree<T>::zMin() const noexcept{ return _root ? _root->zMin() : NAN; }

template <class T>
double Tree<T>::zMax() const noexcept{ return _root ? _root->zMax() : NAN; }

template <class T>
std::vector<std::unique_ptr<Shape>> Tree<T>::remove(T& node){
    PtrList list;
    destruct(node, list);
    return list;
}

template <class T>
std::unique_ptr<Shape> Tree<T>::boundingBox(const_iterator begin, const_iterator end) const noexcept{
    double x0, y0, z0, x1, y1, z1;
    x0 = y0 = z0 = std::numeric_limits<double>::max();
    x1 = y1 = z1 = std::numeric_limits<double>::min();
    for (auto it = begin; it != end; it++){
        x0 = std::min(x0, (*it)->xMin());
        y0 = std::min(y0, (*it)->yMin());
        z0 = std::min(z0, (*it)->zMin());
        x1 = std::max(x1, (*it)->xMax());
        y1 = std::max(y1, (*it)->yMax());
        z1 = std::max(z1, (*it)->zMax());
    }
    const double epsx = eps*(x1-x0);
    const double epsy = eps*(y1-y0);
    const double epsz = eps*(z1-y0);
    x0 -= epsx; y0 -= epsy; z0 -= epsz;
    x1 += epsx; y1 += epsy; z1 += epsz;
    return std::make_unique<Box>(x0, y0, z0, x1, y1, z1);
}

template<typename T>
double Tree<T>::intensity(const Ray& ray, bool isRefracted, bool isReflected) const noexcept{
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
        Shape* next = nextShape(top.position(), top.direction(), top.host(), s);
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
            if (top.host() == next) top.setIntensity(top.intensity()*exp(-next->Sigma_t()*s)); // moving within a Shape, implicit absorption

            if (!isRefracted && !isReflected){
                top.setHost(next);
                stack.push(std::move(top));
            } else{
                double n1, n2;
                if (!top.host()){ // current points to empty space
                    n1 = 1.0;
                    n2 = next->refractive();
                } else if (s == 0.0){
                    // particle enters directly from current to next
                    n1 = top.host()->refractive();
                    n2 = next->refractive();
                } else if (top.host() != next){
                    // particle travels from one Node to another with vacuum in between them
                    n1 = 1.0;
                    n2 = next->refractive();
                } else{
                    // particle travels within the same Node
                    n1 = top.host()->refractive();
                    auto pseudoNext = nextShape(top.position(), top.direction(), top.host(), s);
                    n2 = (s == 0) ? pseudoNext->refractive() : 1.0;
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
template class Tree<Node>;

#include "UNode.h"
template class Tree<UNode>;