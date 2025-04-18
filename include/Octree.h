#ifndef OCTREE_H
#define OCTREE_H
#include "BoundingBox.h"

class Octree{
    public:
    Octree(): _root(nullptr) {}
    explicit Octree(std::vector<Node>& nodes);
    Octree(const Octree&) = delete;
    Octree(Octree&&) = default;
    virtual ~Octree() = default;
    
    Octree& operator=(const Octree&) = delete;
    Octree& operator=(Octree&&) = default;

    void insert(Node& node);
    std::vector<Node> remove(Node& node);
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
    Shape* nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const;

    double xMin() const noexcept;
    double xMax() const noexcept;
    double yMin() const noexcept;
    double yMax() const noexcept;
    double zMin() const noexcept;
    double zMax() const noexcept;

    explicit operator bool() const noexcept{ return bool(_root); } // to check if tree is empty
    friend std::ostream& operator<<(std::ostream& os, const Octree& tree);

    protected:
    Node _root;
};

#endif // OCTREE_H