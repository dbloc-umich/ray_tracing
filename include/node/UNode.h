/**
 * UNode is short for "Node with Unique leaves", which is a simple Node but with an extra set of pointers pointing to its leaves.
 * The leaves for unique children of the box that contains it, not any other Node.
*/

#ifndef UNODE_H
#define UNODE_H

#include "NodeBase.h"
#include <vector>

class Node;
class UNode: public NodeBase<UNode, 8>{
    public:
    using NodeBase<UNode, 8>::NodeBase;

    // Access the leaves
    std::unique_ptr<Node>& operator()(std::size_t i) noexcept{ return _leaves[i]; }
    const std::unique_ptr<Node>& operator()(std::size_t i) const noexcept{ return _leaves[i]; }
    void push(std::unique_ptr<Node> node) noexcept{ _leaves.push_back(std::move(node)); }

    bool leavesOverlap(const Shape& other) const noexcept override;
    bool leavesOverlap(const UNode& other) const noexcept override;

    std::size_t leafCount() const noexcept{ return _leaves.size(); }
    bool empty() const noexcept override;
    std::size_t octant(const Shape& other) const noexcept;

    protected:
    std::vector<std::unique_ptr<Node>> _leaves;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif