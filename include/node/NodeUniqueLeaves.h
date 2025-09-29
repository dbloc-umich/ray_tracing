#ifndef NODE_UNIQUE_LEAVES_H
#define NODE_UNIQUE_LEAVES_H

#include "Node.h"

class NodeUniqueLeaves: public Node{
    public:
    using Node::Node;

    // Access the leaves
    std::unique_ptr<Node>& operator()(std::size_t i) noexcept{ return _leaves[i]; }
    const std::unique_ptr<Node>& operator()(std::size_t i) const noexcept{ return _leaves[i]; }
    void push(std::unique_ptr<Node> node) noexcept{ _leaves.push_back(std::move(node)); }

    bool leavesOverlap(const Shape& other) const noexcept override;
    bool leavesOverlap(const Node& other) const noexcept override;

    std::size_t numContents() const noexcept{ return _leaves.size(); }
    bool empty() const noexcept override;
    std::size_t octant(const Shape& other) const noexcept;

    protected:
    std::vector<std::unique_ptr<Node>> _leaves;
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif