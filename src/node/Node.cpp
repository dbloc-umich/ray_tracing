#include "Node.h"
#include "Shape.h"

Node::Node(std::size_t N):
    _shape(),
    _children(N),
    _level(0)
{}

Node::Node(std::unique_ptr<Shape> shape, std::size_t N):
    _shape(std::move(shape)),
    _children(N),
    _level(0)
{}

bool Node::leavesOverlap(const Shape& other) const noexcept{
    // Check if there's any leaf node that overlaps with the given Shape
    if (!_shape) return false;
    if (!_shape->overlaps(other)) return false;
    for (const auto& it: _children){
        if (!it) continue;
        if (it->isLeaf() && (*it)->overlaps(other)) return true;
        if (it->_shape && it->leavesOverlap(other)) return true;
    }
    return false;
}

bool Node::leavesOverlap(const Node& other) const noexcept{
    if (!_shape || !other._shape) return false;
    if (!_shape->overlaps(*other)) return false;
    for (const auto& it1: _children){
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

std::size_t Node::size() const noexcept{
    std::size_t sz = 0;
    for (const auto& node: _children){
        if (node) sz++;
    }
    return sz;
}

bool Node::empty() const noexcept{
    for (const auto& node: _children){
        if (node) return false;
    }
    return true;
}

std::ostream& Node::printTabs(std::ostream& os, std::size_t n) const noexcept{
    for (unsigned i = 0; i < n; i++) os << "\t";
    return os;
}

std::ostream& Node::print(std::ostream& os) const noexcept{
    printTabs(os, _level);
    os << *_shape << " ";
    
    if (!isLeaf()){
        os << " encloses " << size() << " object";
        if (size() > 1) os << "s";
        os << ":\n";

        for (const auto& it: _children){
            if (it){
                printTabs(os, _level+1);
                os << *it << "\n";
            }
        }
    }
    return os;
}