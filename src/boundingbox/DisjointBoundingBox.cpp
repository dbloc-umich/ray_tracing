// #include "DisjointBoundingBox.h"

// bool DisjointBoundingBox::contentsOverlap(const Shape& other) const noexcept{
//     if (this == &other) return true;  // same object 
//     if (!this->overlaps(other)) return false;

//     for (auto& it: _contents){
//         if (it && it->overlaps(other)) return true;
//     }
//     return BoundingBox::contentsOverlap(other);
// }

// double DisjointBoundingBox::solidVolume() const noexcept{
//     double V = 0.0;
//     for (auto& it: _contents) V += it->volume();
//     for (auto& it: _children){
//         if (it){
//             if (auto box = dynamic_cast<DisjointBoundingBox*>(it.get())) V += box->solidVolume();
//             else V += it->volume();
//         }
//     }
//     return V;
// }

// bool DisjointBoundingBox::empty() const noexcept{
//     for (const auto& node: _contents){
//         if (node) return false;
//     }
//     return BoundingBox::empty();
// }

// std::size_t DisjointBoundingBox::octant(const Shape& other) const noexcept{
//     if (!this->encloses(other)) return -1;
    
//     // Check if there's an octant that fully contains the Shape
//     double xMid = (xMin()+xMax())/2;
//     double yMid = (yMin()+yMax())/2;
//     double zMid = (zMin()+zMax())/2;

//     auto x = other.xMin() < xMid && other.xMax() > xMid;
//     auto y = other.yMin() < yMid && other.yMax() > yMid;
//     auto z = other.zMin() < zMid && other.zMax() > zMid;
//     if (x || y || z) return 8; // should be contained in _contents

//     std::size_t oct = 0;
//     if (other.xMin() >= xMid) oct += 4;
//     if (other.yMin() >= yMid) oct += 2;
//     if (other.zMin() >= zMid) oct += 1;
//     return oct;
// }

// std::ostream& DisjointBoundingBox::print(std::ostream& os) const noexcept{
//     printTabs(os, _level);
//     Box::print(os);
    
//     if (empty()) os << " is empty.";
//     else{
//         os << " encloses " << size() + _contents.size() << " object";
//         if (size() > 1) os << "s";
//         os << ":\n";
//     }

//     for (const auto& it: _contents){
//         if (it){
//             printTabs(os, _level+1);
//             os << *it << "\n";
//         }
//     }
//     for (const auto& it: _children){
//         if (it){
//             if (!dynamic_cast<const BoundingBox*>(it.get())) printTabs(os, _level+1);
//             os << *it << "\n";
//         }
//     }
//     return os;
// }