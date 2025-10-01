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

std::size_t UNode::octant(const Shape& other) const noexcept{
    if (!_shape->encloses(other)) return -1;
    
    // Check if there's an octant that fully contains the Shape
    double xMid = (_shape->xMin()+_shape->xMax())/2;
    double yMid = (_shape->yMin()+_shape->yMax())/2;
    double zMid = (_shape->zMin()+_shape->zMax())/2;

    bool x = other.xMin() < xMid && other.xMax() > xMid;
    bool y = other.yMin() < yMid && other.yMax() > yMid;
    bool z = other.zMin() < zMid && other.zMax() > zMid;
    if (x || y || z) return 8; // should be contained in _leaves

    std::size_t oct = 0;
    if (other.xMin() >= xMid) oct += 4;
    if (other.yMin() >= yMid) oct += 2;
    if (other.zMin() >= zMid) oct += 1;
    return oct;
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