#include "kDTree.h"
#include "Direction.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stack>
#include <tuple>

namespace{
    using iterator = NodeList::iterator;
    using const_iterator = NodeList::const_iterator;
    using BoxStack = std::stack<BoundingBox*>;
    constexpr double eps = 1e-5; // a small nudge given to each dimension of the BoundingBox

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
    x0 -= eps; y0 -= eps; z0 -= eps;
    x1 += eps; y1 += eps; z1 += eps;
    return std::make_tuple(Point(x0, y0, z0), Point(x1, y1, z1));
}

iterator medianNode(BoundingBox& bbox0, BoundingBox& bbox1, iterator begin, iterator end, const char axis='x'){
    /**
     * Inputs:
     *  bbox0 is the bounding box that completely encloses all objects in nodes[begin:end)
     *  bbox1 must be empty
     * Outputs:
     *  since this function is only called when N = end-begin > 8:
     *      returns an iterator pointing to an element N/2 away from begin, where N=end-begin
     *      the underlying nodeList is partially sorted such that the elements in [begin, median) are contained in bbox0, and those in [median, end) are contained in bbox1
     *      bbox0 is now modified to only contain the lower half of nodes0
    **/
    
    auto comp = [axis](const Node& a, const Node& b)
        {
            if (axis == 'x') return a->xMin() < b->xMin();
            else if (axis == 'y') return a->yMin() < b->yMin();
            else return a->zMin() < b->zMin();
        };
    const std::ptrdiff_t N = end-begin;
    iterator median = begin+N/2;
    std::nth_element(begin, median, end, comp); // Find the median
    
    // bbox1 now should encloses the upper half of the nodes in nodes, now partialy sorted about element number N/2-1
    Point p = bbox0.lowerVertex();
    if (axis == 'x') p.setX((*median)->xMin() - eps);
    else if (axis == 'y') p.setY((*median)->yMin() - eps);
    else p.setZ((*median)->zMin() - eps);
    bbox1.setVertices(p, bbox0.upperVertex());

    // Restructure bbox0 to only enclose the lower half of nodes0
    Point lv, uv;
    std::tie(lv, uv) = findVertices(begin, median);
    bbox0.setVertices(lv, uv);
    return median;
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
        if (end - begin <= 8){
            for (auto it = begin; it != end; it++) (*node)[it-begin] = std::move(*it);
        } else{
            std::vector<BoundingBox> bbox(9);
            std::vector<iterator> sep(9, begin);
            bbox[0] = BoundingBox(lower, upper);
            sep[8] = end;

            // Divide a single box along the x-axis
            sep[4] = medianNode(bbox[0], bbox[4], begin, end, 'x');

            // Divide the 2 boxes along the y-axis
            for (std::size_t i = 0; i < 8; i+=4){
                if (sep[i+4]-sep[i] > 8){
                    sep[i+2] = medianNode(bbox[i], bbox[i+2], sep[i], sep[i+4], 'y');
                }
            }

            // Divide the 4 boxes along the z-axis
            for (std::size_t i = 0; i < 8; i+=2){
                bool splitsInZ = (sep[i+2]-sep[i] > (int)8) && (sep[i+2]-sep[i] < end-begin);
                if (splitsInZ) sep[i+1] = medianNode(bbox[i], bbox[i+1], sep[i], sep[i+2], 'z');
            }

            for (std::size_t i = 1; i < 8; ){
                if (sep[i] == begin){
                    bbox.erase(bbox.begin()+i);
                    sep.erase(sep.begin()+i);
                } else i++;
            }
            bbox.pop_back();

            for (std::size_t i = 0; i < bbox.size(); i++){
                Point lv = bbox[i].lowerVertex();
                Point uv = bbox[i].upperVertex();
                construct((*node)[i], sep[i], sep[i+1], lv, uv, level+1);
            } 
        }
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

bool hasValidBoundingBoxes(const Node& current){
    if (!current) return true; // nullptr
    if (const BoundingBox* box = dynamic_cast<const BoundingBox*>(current.get())){
        for (std::size_t i = 0; i < box->size(); i++){
            if (!(*box)[i]) continue; // nullptr
            if (!box->encloses( *(*box)[i] )) return false;
            if (!hasValidBoundingBoxes((*box)[i])) return false;
        }
        return true;
    }
    return true; // not a BoundingBox;
}

bool hasOverlappingContents(const Node& current){
    // Potentially O(N^2) in time complexity, 
    if (!current) return false; // nullptr;
    if (const BoundingBox* box = dynamic_cast<const BoundingBox*>(current.get())){
        for (std::size_t i = 0; i < box->size(); i++){
            if (!(*box)[i]) continue; // nullptr
            if (hasOverlappingContents((*box)[i])) return true;

            for (std::size_t j = i+1; j < box->size(); j++){
                if (!(*box)[j]) continue; // nullptr
                if ((*box)[i]->contentsOverlap(*(*box)[j])) return true;
            }
        }
    }
    return false; // not a BoundingBox
}

Shape* nextExternalNode(BoxStack& stack, const Point& pos, const Direction& dir, Shape* currentNode, double& s){
    if (stack.empty()) return nullptr;

    BoundingBox* box = stack.top();
    std::size_t index = box->max_size();
    s = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < box->size(); i++){
        if ((*box)[i].get() == currentNode) continue;
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
        currentNode = box;
        stack.pop();
        return nextExternalNode(stack, pos, dir, currentNode, s);
    }
    if (BoundingBox* subbox = dynamic_cast<BoundingBox*>((*box)[index].get())){ // Goes deeper into currentNode
        stack.push(subbox);
        return nextExternalNode(stack, pos, dir, currentNode, s);
    }
    // Found it
    return (*box)[index].get();
}
}

kDTree::kDTree(NodeList& nodes):
    _root(nullptr)
{
    if (!nodes.empty()){
        Point lower, upper;
        std::tie(lower, upper) = findVertices(nodes.cbegin(), nodes.cend());
        construct(_root, nodes.begin(), nodes.end(), lower, upper, 0);
        
        // The first throw is to check where construct() is implemented correctly. This throw right here might even be a bad design.
        //if (!hasValidBoundingBoxes(_root)) throw std::runtime_error("ERROR: At least one bounding box does not enclose its contents.");
        if (hasOverlappingContents(_root)) throw std::invalid_argument("ERROR: There are overlapping contents.");
        nodes.clear();
    }
}

void kDTree::insert(Node& node){
    if (_root->contentsOverlap(*node))
        throw std::invalid_argument("ERROR: The node to be inserted overlaps with the content of the tree.");
    // Find the deepest level that node can be inserted
}

NodeList kDTree::remove(Node& node){
    NodeList list;
    destruct(node, list);
    return list;
}

Shape* kDTree::nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const{
    // Check current to see whether the Point is inside the Shape or it is entering from the surface
    if (current){
        s = current->distanceToSurface(pos, dir);
        if (current->encloses(pos) || (current->surfaceContains(pos) && s > 0.0)) return current;
    }

    // current is guaranteed to not be the returned pointer
    if (_root && !std::isnan(_root->distanceToSurface(pos, dir))){
        BoxStack stack;
        stack.emplace(dynamic_cast<BoundingBox*>(_root.get()));
        auto nextNode = nextExternalNode(stack, pos, dir, current, s);
        if (!nextNode) s = NAN;
        return nextNode;
    }
    s = NAN;
    return nullptr;
}

double kDTree::xMin() const noexcept{ return _root ? _root->xMin() : NAN; }
double kDTree::xMax() const noexcept{ return _root ? _root->xMax() : NAN; }
double kDTree::yMin() const noexcept{ return _root ? _root->yMin() : NAN; }
double kDTree::yMax() const noexcept{ return _root ? _root->yMax() : NAN; }
double kDTree::zMin() const noexcept{ return _root ? _root->zMin() : NAN; }
double kDTree::zMax() const noexcept{ return _root ? _root->zMax() : NAN; }

std::ostream& operator<<(std::ostream& os, const kDTree& tree){
    if (tree._root) os << *(tree._root);
    else os << "Empty tree.";
    return os;
}