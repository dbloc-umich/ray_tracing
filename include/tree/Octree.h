// #ifndef OCTREE_H
// #define OCTREE_H

// #include "Tree.h"
// #include "UNode.h"

// class Octree: public Tree<UNode>{
//     public:
//     using Tree<UNode>::Tree;
//     explicit Octree(PtrList& ptrs);

//     void insert(std::unique_ptr<Shape> shape) override;
//     Shape* nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const noexcept override;

//     protected:
//     void construct(UNode& current, PtrList& ptrs, std::size_t level=0) noexcept;
//     void destruct(UNode& current, PtrList& ptrs) noexcept override;
//     UNode& smallestBox(UNode& current, const std::unique_ptr<Shape>& shape) const noexcept override;
//     bool hasOverlappingContents(const UNode& current) const noexcept override;
// };

// #endif // OCTREE_H