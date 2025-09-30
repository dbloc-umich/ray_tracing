/**
 * This is a simple Node class, containing a pointer to a Shape and an array of 2 pointers to its child Nodes.
*/

#ifndef NODE_H
#define NODE_H

#include "NodeBase.h"

class Shape;
class Node: public NodeBase<Node, 2>{
    public:
    using NodeBase<Node, 2>::NodeBase;

    bool leavesOverlap(const Shape& other) const noexcept override;
    bool leavesOverlap(const Node& other) const noexcept override;
    bool empty() const noexcept override;

    protected:
    std::ostream& print(std::ostream& os) const noexcept override;
};

#endif