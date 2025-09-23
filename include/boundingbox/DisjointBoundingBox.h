// /**
//  * A disjoint bounding box is one where none of its children overlaps.
//  * Because of this property, its eight children are separated exactly at the midpoint of each axis.
//  * Any leaf node that crosses a midpoint on any axis is contained within a separate _contents vector.
// */

// #ifndef DBOUNDING_BOX_H
// #define DBOUNDING_BOX_H
// #include "BoundingBox.h"
// #include <vector>

// class DisjointBoundingBox: public BoundingBox<8>{
//     public:
//     using BoundingBox::BoundingBox;

//     // Access the contents
//     Node& operator()(std::size_t i) noexcept{ return _contents[i]; }
//     const Node& operator()(std::size_t i) const noexcept{ return _contents[i]; }
//     void push(Node& node) noexcept{ _contents.push_back(std::move(node)); }
//     void push(Node&& node) noexcept{ _contents.push_back(std::move(node)); }

//     bool contentsOverlap(const Shape& other) const noexcept override;

//     double solidVolume() const noexcept override;
//     std::size_t numContents() const noexcept{ return _contents.size(); }
//     bool empty() const noexcept override;
//     std::size_t octant(const Shape& other) const noexcept;

//     private:
//     std::vector<Node> _contents{};
//     std::ostream& print(std::ostream& os) const noexcept override;
// };
// #endif // DBOUNDING_BOX_H