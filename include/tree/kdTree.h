#ifndef KDTREE_H
#define KDTREE_H
#include "Tree.h"

class kdTree: public Tree{
    public:
    using Tree::Tree;
    explicit kdTree(std::vector<Node>& nodes);

    void insert(Node& node) override;
    Shape* nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const override;

    protected:
    iterator medianNode(BoundingBox& bbox0, BoundingBox& bbox1, iterator begin, iterator end, const char axis);
    void construct(Node& current, iterator begin, iterator end, Point& lower, Point& upper, std::size_t level);
    void destruct(Node& current, NodeList& nodes) override;
    bool hasOverlappingContents(const Node& current) const override;
};

#endif // KDTREE_H