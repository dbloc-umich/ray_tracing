#ifndef KD_TREE_H
#define KD_TREE_H
#include "Tree.h"
class kDTree : public Tree{
    public:
    kDTree(): Tree() {}
    explicit kDTree(NodeList& nodes);

    void insert(Node& node) override;
    NodeList remove(Node& node) override;
    Shape* nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const override;
};
#endif // KD_TREE_H