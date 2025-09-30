#include "UNode.h"
#include "Node.h"
#include "Shape.h"

bool UNode::leavesOverlap(const Shape& other) const noexcept{
    // Check if there's any leaf node that overlaps with the given Shape
    if (!_shape) return false;
    if (!_shape->overlaps(other)) return false;
    for (auto& it: _leaves){
        if ((*it)->overlaps(other)) return true;
    }
    for (auto& it: _children){
        if (!it) continue;
        if (it->isLeaf() && (*it)->overlaps(other)) return true;
        if (it->_shape && it->leavesOverlap(other)) return true;
    }
    return false;
}

bool UNode::leavesOverlap(const UNode& other) const noexcept{
    if (!_shape || !other._shape) return false;
    if (!_shape->overlaps(*other)) return false;
    for (auto& it: _leaves){
        if (other.leavesOverlap(**it)) return true;
    }
    for (auto& it1: _children){
        if (!it1) continue;
        bool isLeaf1 = it1->isLeaf();
        for (const auto& it2: other._children){
            // it is std::unique_ptr<Node>, *it is Node&, and **it is Shape&
            if (!it2) continue;
            bool isLeaf2 = it2->isLeaf();
            if (isLeaf1 && isLeaf2){
                if ((*it1)->overlaps(**it2)) return true;
            } else if (isLeaf1){
                if (it2->leavesOverlap(**it1)) return true;
            } else if (isLeaf2){
                if (it1->leavesOverlap(**it2)) return true;
            } else if (it1->leavesOverlap(*it2)) return true;
        }
    }
    return false;
}

bool UNode::empty() const noexcept{
    bool x = std::find_if(_children.cbegin(), _children.cend(), [](auto& node){ return bool(node); }) == _children.cend();
    bool y = _leaves.empty();
    return x && y;
}

std::ostream& UNode::print(std::ostream& os) const noexcept{
    printTabs(os, _level);
    os << *_shape << " ";
    
    if (!isLeaf()){
        os << " encloses " << size()+leafCount() << " object";
        if (size()+leafCount() > 1) os << "s";
        os << ":\n";

        for (auto& it: _leaves){
            if (it){
                printTabs(os, _level+1);
                os << *it << "\n";
            }
        }

        for (auto& it: _children){
            if (it){
                printTabs(os, _level+1);
                os << *it << "\n";
            }
        }
    }
    return os;
}