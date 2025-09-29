// #ifndef OCTREE_H
// #define OCTREE_H
// #include "Tree.h"

// class Octree: public Tree{
//     public:
//     using Tree::Tree;
//     explicit Octree(std::vector<Node>& nodes);

//     void insert(Node& node) override {};
//     Shape* nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const override;

//     protected:
//     void construct(Node& current, PtrList& nodes, const Point& lower, const Point& upper, std::size_t level=0);
//     void destruct(Node& current, PtrList& nodes) override;
//     bool hasOverlappingContents(const Node& current) const override;
// };

// #endif // OCTREE_H