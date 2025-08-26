#include "BoundingBox.h"
#include <string>

static std::ostream& printTabs(std::ostream& os, std::size_t n){
    for (unsigned i = 0; i < n; i++) os << "\t";
    return os;
}

BoundingBox::BoundingBox(const Point& lower, const Point& upper):
    Box(lower, upper),
    _children{},
    _contents{},
    _level(0)
{}

BoundingBox::BoundingBox(double x0, double y0, double z0, double x1, double y1, double z1):
    Box(x0, y0, z0, x1, y1, z1),
    _children{},
    _contents{},
    _level(0)    
{}

bool BoundingBox::contentsOverlap(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object 
    if (!this->overlaps(other)) return false;

    for (auto& it: _children){
        if (it && it->contentsOverlap(other)) return true;
    }
    for (auto& it: _contents){
        if (it && it->overlaps(other)) return true;
    }
    return false;
}

std::size_t BoundingBox::size() const noexcept{
    std::size_t sz = 0;
    for (const auto& node: _children){
        if (node) sz++;
    }
    return sz;
}

bool BoundingBox::empty() const noexcept{
    for (const auto& node: _children){
        if (node) return false;
    }
    for (const auto& node: _contents){
        if (node) return false;
    }
    return true;
}

bool BoundingBox::full() const noexcept{
    for (const auto& node: _children){
        if (!node) return false;
    }
    return true;    
}

std::size_t BoundingBox::octant(const Shape& other) const noexcept{
    if (!this->encloses(other)) return -1;
    
    // Check if there's an octant that fully contains the Shape
    double xMid = (xMin()+xMax())/2;
    double yMid = (yMin()+yMax())/2;
    double zMid = (zMin()+zMax())/2;

    auto x = other.xMin() < xMid && other.xMax() > xMid;
    auto y = other.yMin() < yMid && other.yMax() > yMid;
    auto z = other.zMin() < zMid && other.zMax() > zMid;
    if (x || y || z) return 8; // should be contained in _contents

    std::size_t oct = 0;
    if (other.xMin() >= xMid) oct += 4;
    if (other.yMin() >= yMid) oct += 2;
    if (other.zMin() >= zMid) oct += 1;
    return oct;
}

std::ostream& BoundingBox::print(std::ostream& os) const noexcept{
    printTabs(os, _level);
    Box::print(os);
    
    if (empty()) os << " is empty.";
    else{
        os << " encloses " << size() + _contents.size() << " object";
        if (size() > 1) os << "s";
        os << ":\n";
    }

    for (const auto& it: _contents){
        if (it){
            printTabs(os, _level+1);
            os << *it << "\n";
        }
    }
    for (const auto& it: _children){
        if (it){
            if (!dynamic_cast<const BoundingBox*>(it.get())) printTabs(os, _level+1);
            os << *it << "\n";
        }
    }
    return os;
}