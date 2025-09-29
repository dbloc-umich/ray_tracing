#ifndef TREE_H
#define TREE_H

#include "Node.h"

class Box;
class Direction;
class Point;
class Tree{
    public:
    Tree(): _root(nullptr) {}
    Tree(const Tree&) = delete;
    Tree(Tree&&) = default;
    virtual ~Tree() = default;
    Tree& operator=(const Tree&) = delete;
    Tree& operator=(Tree&&) = default;

    virtual void insert(std::unique_ptr<Shape> shape) = 0;
    std::vector<std::unique_ptr<Shape>> remove(Node& node);

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

    double xMin() const noexcept;
    double xMax() const noexcept;
    double yMin() const noexcept;
    double yMax() const noexcept;
    double zMin() const noexcept;
    double zMax() const noexcept;

    explicit operator bool() const noexcept{ return bool(_root); } // to check if tree is empty
    Node& root() noexcept{ return _root; }
    const Node& root() const noexcept{ return _root; }
    friend std::ostream& operator<<(std::ostream& os, const Tree& tree);

    protected:
    Node _root;

    static constexpr double eps = 1e-6;
    using PtrList = std::vector<std::unique_ptr<Shape>>;
    using iterator = PtrList::iterator;
    using const_iterator = PtrList::const_iterator;
    std::unique_ptr<Box> boundingBox(const_iterator begin, const_iterator end) const;

    // Required helper functions
    virtual void destruct(Node& current, PtrList& nodes) = 0;
    virtual Node& smallestBox(Node& current, const std::unique_ptr<Shape>& shape) const noexcept = 0;
    virtual bool hasOverlappingContents(const Node& current) const = 0;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Tree& tree);

#endif // TREE_H