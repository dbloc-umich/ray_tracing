#include "BoundingBox.h"
#include <string>

BoundingBox::BoundingBox(const Point& lower, const Point& upper):
    Box(lower, upper),
    _children{},
    _level(0)
{}

BoundingBox::BoundingBox(double x0, double y0, double z0, double x1, double y1, double z1):
    Box(x0, y0, z0, x1, y1, z1),
    _children{},
    _level(0)    
{}

bool BoundingBox::contentsOverlap(const Shape& other) const noexcept{
    if (this == &other) return true;  // same object 
    if (!this->overlaps(other)) return false;

    for (auto& it: _children){
        if (it && it->contentsOverlap(other)) return true;
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

std::ostream& BoundingBox::printTabs(std::ostream& os, std::size_t n) const noexcept{
    for (unsigned i = 0; i < n; i++) os << "\t";
    return os;
}

std::ostream& BoundingBox::print(std::ostream& os) const noexcept{
    printTabs(os, _level);
    Box::print(os);
    
    if (empty()) os << " is empty.";
    else{
        os << " encloses " << size() << " object";
        if (size() > 1) os << "s";
        os << ":\n";
    }

    for (const auto& it: _children){
        if (it){
            if (!dynamic_cast<const BoundingBox*>(it.get())) printTabs(os, _level+1);
            os << *it << "\n";
        }
    }
    return os;
}