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
    _children{},
    _level(0)
{}

BoundingBox::BoundingBox(double x0, double y0, double z0, double x1, double y1, double z1):
    Box(x0, y0, z0, x1, y1, z1),
    _children{},
    _level(0)
{}

BoundingBox::BoundingBox(NodeList& values):
    Box(),
    _children{},
    _level(0)
{
    if (values.size() > 8){
        std::string msg = "ERROR: Each bounding box can only contain at most " + std::to_string(8) + " children.";
        throw std::out_of_range(msg);
    }
    if (values.size() > 0){
        double x0, y0, z0, x1, y1, z1;
        x0 = y0 = z0 = std::numeric_limits<double>::max();
        x1 = y1 = z1 = std::numeric_limits<double>::min();
        for (auto it = values.begin(); it != values.end(); it++){
            if (*it){ // *it is not nullptr
                x0 = std::min(x0, (*it)->xMin());
                y0 = std::min(y0, (*it)->yMin());
                z0 = std::min(z0, (*it)->zMin());
                x1 = std::max(x1, (*it)->xMax());
                y1 = std::max(y1, (*it)->yMax());
                z1 = std::max(z1, (*it)->zMax());
                //(*it)->setParent(this);
                _children.emplace_back(std::move(*it));
            }
        }
        setVertices(Point(x0, y0, z0), Point(x1, y1, z1));
    }
}

Node& BoundingBox::operator[](std::size_t i) noexcept {
    if (i >= _children.size() && i < 8) _children.resize(i+1);
    return _children[i]; // undefined behavior, but consistent with how operator[] is defined in std::vector
}

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

std::ostream& BoundingBox::print(std::ostream& os) const noexcept{
    printTabs(os, _level);
    Box::print(os);
    //std::size_t sz = size();
    
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