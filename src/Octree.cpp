#include "Octree.h"
#include <cmath>

namespace{
    using iterator = NodeList::iterator;
    using const_iterator = NodeList::const_iterator;
    using BoxStack = std::stack<BoundingBox*>;

decltype(auto) findVertices(const_iterator cbegin, const_iterator cend){
    double x0, y0, z0, x1, y1, z1;
    x0 = y0 = z0 = std::numeric_limits<double>::max();
    x1 = y1 = z1 = std::numeric_limits<double>::min();
    for (auto it = cbegin; it != cend; it++){
        x0 = std::min(x0, (*it)->xMin());
        y0 = std::min(y0, (*it)->yMin());
        z0 = std::min(z0, (*it)->zMin());
        x1 = std::max(x1, (*it)->xMax());
        y1 = std::max(y1, (*it)->yMax());
        z1 = std::max(z1, (*it)->zMax());
    }
    x0 -= BOX_EPS; y0 -= BOX_EPS; z0 -= BOX_EPS;
    x1 += BOX_EPS; y1 += BOX_EPS; z1 += BOX_EPS;
    return std::make_tuple(Point(x0, y0, z0), Point(x1, y1, z1));
}

void construct(Node& current, iterator begin, iterator end, Point& lower, Point& upper, std::size_t level){
    /**
     * Inputs:
     *  current: pointer to the BoundingBox at which space partitioning occurs
     *  nodes: a vector of all nodes contained within *current
     *  begin, end: indices of nodes to partition
     *  lower, upper: vertices of the BoundingBox that must enclose all objects in nodes
     * Outputs:
     *  Recursively construct a kDTree as an octree
    **/

    if (!current){ // current is nullptr
        std::unique_ptr<BoundingBox> node = std::make_unique<BoundingBox>(lower, upper);
        node->setLevel(level);
        current = std::move(node);
    }    
}

void destruct(Node& current, NodeList& nodes){
    if (current){ // not nullptr
        if (BoundingBox* box = dynamic_cast<BoundingBox*>(current.get())){
            // If current is a BoundingBox, recursively destruct its children
            for (std::size_t i = 0; i < box->size(); i++) destruct((*box)[i], nodes);
        } else nodes.emplace_back(std::move(current));
        current.reset();
    }
}

bool hasOverlappingContents(const Node& current){
    // Potentially O(N^2) in time complexity

    if (!current) return false; // nullptr;
    if (const BoundingBox* box = dynamic_cast<const BoundingBox*>(current.get())){
        for (std::size_t i = 0; i < box->size(); i++){
            if (!(*box)[i]) continue; // nullptr
            if (hasOverlappingContents((*box)[i])) return true;

            for (std::size_t j = 0; j < box->numContents(); j++){
                if ((*box)[i]->contentsOverlap(*(*box)(j))) return true;
            }
        }
    }
    return false; // not a BoundingBox
}

Shape* nextExternalNode(BoxStack& stack, const Point& pos, const Direction& dir, std::vector<Shape*>& visitedNodes, double& s){ return visitedNodes[0]; }
}

Shape* Octree::nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const{
    // Check current to see whether the Point is inside the Shape or it is entering from the surface
    if (current){
        s = current->distanceToSurface(pos, dir);
        if (current->encloses(pos) || (current->surfaceContains(pos) && s > 0.0)) return current;
    }

    // current is guaranteed to not be the returned pointer
    if (_root && !std::isnan(_root->distanceToSurface(pos, dir))){
        BoxStack stack;
        stack.emplace(dynamic_cast<BoundingBox*>(_root.get()));
        std::vector<Shape*> visited{current};
        visited.reserve(8);
        auto nextNode = nextExternalNode(stack, pos, dir, visited, s);
        if (!nextNode) s = NAN;
        return nextNode;
    }
    s = NAN;
    return nullptr;
}

double Octree::xMin() const noexcept{ return _root ? _root->xMin() : NAN; }
double Octree::xMax() const noexcept{ return _root ? _root->xMax() : NAN; }
double Octree::yMin() const noexcept{ return _root ? _root->yMin() : NAN; }
double Octree::yMax() const noexcept{ return _root ? _root->yMax() : NAN; }
double Octree::zMin() const noexcept{ return _root ? _root->zMin() : NAN; }
double Octree::zMax() const noexcept{ return _root ? _root->zMax() : NAN; }

std::ostream& operator<<(std::ostream& os, const Octree& tree){
    if (tree._root) os << *(tree._root);
    else os << "Empty tree.";
    return os;
}