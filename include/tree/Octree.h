#ifndef OCTREE_H
#define OCTREE_H
#include "Tree.h"

class Octree: public Tree{
    public:
    using Tree::Tree;
    explicit Octree(std::vector<Node>& nodes);

    void insert(Node& node) override;
    Shape* nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const override;

    protected:
    void construct(Node& current, NodeList& nodes, const Point& lower, const Point& upper, std::size_t level=0);
    void destruct(Node& current, NodeList& nodes) override;
    Node& smallestParentNode(const Node& node, Node& parent, Node& current);
    bool hasOverlappingContents(const Node& current) const override;
};

#endif // OCTREE_H