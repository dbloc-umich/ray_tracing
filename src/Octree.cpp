#include "Octree.h"
#include "Direction.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <stack>

namespace{
    static constexpr double eps = 1e-5; // a small nudge given to each dimension of the BoundingBox
    using NodeList = std::vector<Node>;
    using iterator = NodeList::iterator;
    using const_iterator = NodeList::const_iterator;
    using BoxStack = std::stack<BoundingBox*>;

void construct(Node& current, NodeList& nodes, const Point& lower, const Point& upper, std::size_t level=0){

    if (!current){ // current is nullptr
        std::unique_ptr<BoundingBox> box = std::make_unique<BoundingBox>(lower, upper);
        box->setLevel(level);
        if (nodes.size() <= 8) {
            for (std::size_t i = 0; i < nodes.size(); i++) (*box)[i] = std::move(nodes[i]);
        }
        else{
            std::array<NodeList, 8> octantList;
            while (!nodes.empty()){
                const std::size_t oct = box->octant(*(nodes.back()));
                if (oct == 8) box->push(nodes.back());
                else octantList[oct].push_back(std::move(nodes.back()));
                nodes.pop_back();              
            }

            for (std::size_t oct = 0; oct < 8; oct++){
                Point subLower(lower);
                Point subUpper(upper);
                double xMid = (lower.x()+upper.x())/2;
                double yMid = (lower.y()+upper.y())/2;
                double zMid = (lower.z()+upper.z())/2;

                if (oct < 4) subUpper.setX(xMid);
                else subLower.setX(xMid);
                if (oct % 4 < 2) subUpper.setY(yMid);
                else subLower.setY(yMid);
                if (!(oct % 2)) subUpper.setZ(zMid);
                else subLower.setZ(zMid);
                construct((*box)[oct], octantList[oct], subLower, subUpper, level+1);
            }
        }
        current = std::move(box);
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
        for (std::size_t j = 0; j < box->numContents(); j++){
            for (std::size_t k = j+1; k < box->numContents(); k++){
                if ((*box)(j)->overlaps(*(*box)(k))) return true;
            }
        }
        
        for (std::size_t i = 0; i < box->size(); i++){
            for (std::size_t j = 0; j < box->numContents(); j++){
                if ((*box)[i]->contentsOverlap(*(*box)(j))) return true;
            }
            if (hasOverlappingContents((*box)[i])) return true;
        }
    }
    return false; // not a BoundingBox
}

Shape* nextExternalNode(BoxStack& stack, const Point& pos, const Direction& dir, std::vector<Shape*>& visitedNodes, double& s){
    if (stack.empty()) return nullptr;

    BoundingBox* box = stack.top(); // stack.top() is always at least the parent of all nodes in visitedNodes
    std::size_t index = box->max_size();
    s = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < box->size(); i++){
        if (std::find(visitedNodes.cbegin(), visitedNodes.cend(), (*box)[i].get()) != visitedNodes.cend()) continue; // currently on a visited Node, will skip
        double dist = (*box)[i]->distanceToSurface(pos, dir);
        if (!std::isnan(dist)){
            Point nextPos = pos;
            nextPos.advance(dir, dist); // to check against tangency
            if (!dir.isOrthogonal((*box)[i]->normal(nextPos)) && dist < s){
                index = i;
                s = dist;
            }
        }
    }
    
    if (index == box->max_size()){ // Can't find nearest neighbor at currentNode
        if (visitedNodes.size() == stack.top()->size()){
            visitedNodes.clear();
            visitedNodes.push_back(stack.top());
            stack.pop();
        }
        else if (visitedNodes.size() > 1) visitedNodes.push_back(box);
        else{
            if (!dynamic_cast<BoundingBox*>(visitedNodes[0])){
                stack.pop();
                visitedNodes[0] = box;
            }
            else visitedNodes.push_back(box);
        }
        return nextExternalNode(stack, pos, dir, visitedNodes, s);
    }
    if (BoundingBox* subbox = dynamic_cast<BoundingBox*>((*box)[index].get())){ // Goes deeper into currentNode
        stack.push(subbox);
        return nextExternalNode(stack, pos, dir, visitedNodes, s);
    }
    return (*box)[index].get(); // Found it
}
}

Octree::Octree(NodeList& nodes):
    _root(nullptr)
{
    if (!nodes.empty()){
        double x0, y0, z0, x1, y1, z1;
        x0 = y0 = z0 = std::numeric_limits<double>::max();
        x1 = y1 = z1 = std::numeric_limits<double>::min();
        for (const auto& it: nodes){
            x0 = std::min(x0, it->xMin());
            y0 = std::min(y0, it->yMin());
            z0 = std::min(z0, it->zMin());
            x1 = std::max(x1, it->xMax());
            y1 = std::max(y1, it->yMax());
            z1 = std::max(z1, it->zMax());
        }
        const double epsx = eps*(x1-x0);
        const double epsy = eps*(y1-y0);
        const double epsz = eps*(z1-y0);
        x0 -= epsx; y0 -= epsy; z0 -= epsz;
        x1 += epsx; y1 += epsy; z1 += epsz;
        
        Point lower(x0, y0, z0);
        Point upper(x1, y1, z1);
        construct(_root, nodes, lower, upper);
        
        if (hasOverlappingContents(_root)){
            NodeList movedNodes;
            destruct(_root, movedNodes);
            while (!movedNodes.empty()){
                nodes.push_back(std::move(movedNodes.back()));
                movedNodes.pop_back();
            }
            throw std::invalid_argument("ERROR: There are overlapping contents.");
        }
        nodes.clear();
    }
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