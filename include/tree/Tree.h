#ifndef TREE_H
#define TREE_H
#include <vector>
#include "BoundingBox.h"

class Tree{
    public:
    Tree(): _root(nullptr) {}
    Tree(const Tree&) = delete;
    Tree(Tree&&) = default;
    virtual ~Tree() = default;
    Tree& operator=(const Tree&) = delete;
    Tree& operator=(Tree&&) = default;

    virtual void insert(Node& node) = 0;
    std::vector<Node> remove(Node& node);

    /**
     * Inputs:
     *  pos: the position of the particle
     *  dir: the direction of the particle
     *  current: a raw pointer to the Shape that has the particle on its surface
     *  s: a placeholder double&
     * Outputs:
     *  s is the modified distance that the Point has to travel to reach the next Shape
     *  a raw pointer to Shape (from Node::get()) where the particle lands is returned
    **/
    virtual Shape* nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const = 0;

    double xMin() const noexcept;
    double xMax() const noexcept;
    double yMin() const noexcept;
    double yMax() const noexcept;
    double zMin() const noexcept;
    double zMax() const noexcept;

    explicit operator bool() const noexcept{ return bool(_root); } // to check if tree is empty
    Shape* root() const noexcept{ return _root.get(); }
    friend std::ostream& operator<<(std::ostream& os, const Tree& tree);

    protected:
    Node _root;

    static constexpr double eps = 1e-6;
    using NodeList = std::vector<Node>;
    using iterator = NodeList::iterator;
    using const_iterator = NodeList::const_iterator;
    std::pair<Point, Point> findVertices(const_iterator begin, const_iterator end) const;

    // Required helper functions
    virtual void destruct(Node& current, NodeList& nodes) = 0;
    virtual bool hasOverlappingContents(const Node& current) const = 0;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Tree& tree);

#endif // TREE_H