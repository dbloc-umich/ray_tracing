#ifndef KDTREE_H
#define KDTREE_H

#include "Node.h"
#include "Tree.h"

class Box;
class kdTree: public Tree<Node>{
    public:
    using Tree<Node>::Tree;
    explicit kdTree(PtrList& ptrs);

    void insert(std::unique_ptr<Shape> shape) override;
    Shape* nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const override;

    protected:
    iterator medianNode(Box& bbox0, Box& bbox1, iterator begin, iterator end, char axis);
    void construct(Node& current, iterator begin, iterator end, std::size_t level, char axis);
    void destruct(Node& current, PtrList& ptrs) override;
    Node& smallestBox(Node& current, const std::unique_ptr<Shape>& shape) const noexcept override;
    bool hasOverlappingContents(const Node& current) const override;
};

#endif // KDTREE_H