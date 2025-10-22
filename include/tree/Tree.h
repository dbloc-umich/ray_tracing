#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <memory>
#include <vector>

class Direction;
class Point;
class Ray;
class Shape;

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
    virtual Shape* nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const noexcept = 0;

    double xMin() const noexcept;
    double xMax() const noexcept;
    double yMin() const noexcept;
    double yMax() const noexcept;
    double zMin() const noexcept;
    double zMax() const noexcept;

    T& root() noexcept{ return _root; }
    const T& root() const noexcept{ return _root; }
    explicit operator bool() const noexcept{ return bool(_root); } // to check if tree is empty
    bool leavesOverlap(const Shape& other) const noexcept{ return _root.leavesOverlap(other); }
    bool leavesOverlap(const T& other) const noexcept{ return _root.leavesOverlap(other); }
   
    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Tree<U>& tree);

    protected:
    T _root;

    static constexpr double eps = 1e-6;
    using PtrList = std::vector<std::unique_ptr<Shape>>;
    using iterator = PtrList::iterator;
    using const_iterator = PtrList::const_iterator;
    std::unique_ptr<Shape> boundingBox(const_iterator begin, const_iterator end) const noexcept;

    // Required helper functions
    virtual void destruct(T& current, PtrList& nodes) noexcept = 0;
    virtual T& smallestBox(T& current, const std::unique_ptr<Shape>& shape) const noexcept = 0;
    virtual bool hasOverlappingContents(const T& current) const noexcept = 0;
};

template <typename U>
std::ostream& operator<<(std::ostream& os, const Tree<U>& tree){
    if (tree._root) os << tree._root;
    else os << "Empty tree.";
    return os;
}

#endif // TREE_H