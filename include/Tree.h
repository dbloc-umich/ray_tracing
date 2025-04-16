#ifndef TREE_H
#define TREE_H
#include "BoundingBox.h"
#include <stack>

static constexpr double BOX_EPS = 1e-5; // a small nudge given to each dimension of the BoundingBox
class Tree{
    public:
    Tree(): _root(nullptr) {}
    Tree(const Tree&) = delete;
    Tree(Tree&&) = default;
    virtual ~Tree() = default;
    
    Tree& operator=(const Tree&) = delete;
    Tree& operator=(Tree&&) = default;

    virtual void insert(Node& node) = 0;
    virtual NodeList remove(Node& node) = 0;
    /**
     * Inputs:
     *  pos: the position of the particle, must be on the surface of a leaf node
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
    friend std::ostream& operator<<(std::ostream& os, const Tree& tree);

    using iterator = NodeList::iterator;
    using const_iterator = NodeList::const_iterator;
    using BoxStack = std::stack<BoundingBox*>;

    protected:
    Node _root;
};

inline decltype(auto) findVertices(Tree::const_iterator cbegin, Tree::const_iterator cend);

#endif // TREE_H