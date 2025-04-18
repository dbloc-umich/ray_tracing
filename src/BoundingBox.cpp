#include "BoundingBox.h"
#include <exception>
#include <limits>
#include <string>

static std::ostream& printTabs(std::ostream& os, std::size_t n){
    for (unsigned i = 0; i < n; i++) os << "\t";
    return os;
}

BoundingBox::BoundingBox(const Point& lower, const Point& upper):
    Box(lower, upper),
    _children(8),
    _contents{},
    _level(0)
{}

BoundingBox::BoundingBox(double x0, double y0, double z0, double x1, double y1, double z1):
    Box(x0, y0, z0, x1, y1, z1),
    _children(8),
    _contents{},
    _level(0)    
{}

bool BoundingBox::contentsOverlap(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object 
    if (!this->overlaps(other)) return false;

    if (const BoundingBox* box = dynamic_cast<const BoundingBox*>(&other)){ // other is also a BoundingBox
        for (auto& it1: _children){
            if (!it1) continue;
            if (dynamic_cast<const BoundingBox*>(it1.get())){ // it1 points to a BoundingBox
                for (auto& it2: box->_children){
                    if (it2 && it1->contentsOverlap(*it2)) return true;
                }
            } else{ // it1 doesn't point to a BoundingBox
                for (auto& it2: box->_children){
                    if (it2 && it2->contentsOverlap(*it1)) return true;
                }
            }
        }
    } else{ // other is not a BoundingBox
        for (auto& it: _children){
            if (it && it->contentsOverlap(other)) return true;
        }
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
    return true;
}

bool BoundingBox::full() const noexcept{
    for (const auto& node: _children){
        if (!node) return false;
    }
    return true;    
}

std::size_t BoundingBox::octant(const Shape& other) const noexcept{
    // Check if there's an octant that fully contains the Shape
    
    std::size_t oct = 0;
    return oct;
}

std::ostream& BoundingBox::print(std::ostream& os) const noexcept{
    printTabs(os, _level);
    Box::print(os);
    
    if (size() == 0) os << " is empty.";
    else{
        os << " encloses " << size() << " object";
        if (size() > 1) os << "s";
        os << ":\n";
    }

    for (auto it = _children.cbegin(); it != _children.cend(); it++){
        if (*it){
            if (!dynamic_cast<const BoundingBox*>(it->get())) printTabs(os, _level+1);
            os << **it << "\n";
        }
    }
    return os;
}