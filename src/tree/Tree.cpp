#include "Tree.h"

#include <cmath>
#include <limits>

#include "Box.h"

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

// Explicit instantiations
#include "Node.h"
template class Tree<Node>;

#include "UNode.h"
template class Tree<UNode>;