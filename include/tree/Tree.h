#ifndef TREE_H
#define TREE_H

#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include "Point.h"
#include "Shape.h"

class Box;
class Direction;
class Point;

template <class T>
class Tree{
    public:
    Tree(): _root(nullptr) {}
    Tree(const Tree<T>&) = delete;
    Tree(Tree<T>&&) = default;
    virtual ~Tree() = default;
    Tree<T>& operator=(const Tree<T>&) = delete;
    Tree<T>& operator=(Tree<T>&&) = default;

    virtual void insert(std::unique_ptr<Shape> shape) = 0;
    std::vector<std::unique_ptr<Shape>> remove(T& node);

    /**
     * Inputs:
     *  pos: the position of the particle
     *  dir: the direction of the particle
     *  current: a raw pointer to the Shape that has the particle on its surface
     *  s: a placeholder double&
     * Outputs:
     *  s is the modified distance that the Point has to travel to reach the next Shape
     *  a raw pointer to Shape where the particle lands is returned
    **/
    virtual Shape* nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const = 0;

    double xMin() const noexcept{ return _root ? _root->xMin() : NAN; }
    double xMax() const noexcept{ return _root ? _root->xMax() : NAN; }
    double yMin() const noexcept{ return _root ? _root->yMin() : NAN; }
    double yMax() const noexcept{ return _root ? _root->yMax() : NAN; }
    double zMin() const noexcept{ return _root ? _root->zMin() : NAN; }
    double zMax() const noexcept{ return _root ? _root->zMax() : NAN; }

    explicit operator bool() const noexcept{ return bool(_root); } // to check if tree is empty
    T& root() noexcept{ return _root; }
    const T& root() const noexcept{ return _root; }
    
    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Tree<U>& tree);

    protected:
    T _root;

    static constexpr double eps = 1e-6;
    using PtrList = std::vector<std::unique_ptr<Shape>>;
    using iterator = PtrList::iterator;
    using const_iterator = PtrList::const_iterator;
    std::unique_ptr<Box> boundingBox(const_iterator begin, const_iterator end) const;

    // Required helper functions
    virtual void destruct(T& current, PtrList& nodes) = 0;
    virtual T& smallestBox(T& current, const std::unique_ptr<Shape>& shape) const noexcept = 0;
    virtual bool hasOverlappingContents(const T& current) const = 0;
};

template <class T>
std::vector<std::unique_ptr<Shape>> Tree<T>::remove(T& node){
    PtrList list;
    destruct(node, list);
    return list;
}

template <typename U>
std::ostream& operator<<(std::ostream& os, const Tree<U>& tree){
    if (tree._root) os << tree._root;
    else os << "Empty tree.";
    return os;
}

template <class T>
std::unique_ptr<Box> Tree<T>::boundingBox(const_iterator begin, const_iterator end) const{
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

#endif // TREE_H