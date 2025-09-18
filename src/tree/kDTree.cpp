#include "kdTree.h"
#include "Direction.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <stack>

namespace{
    using BoxType = BoundingBox;
}

kdTree::kdTree(NodeList& nodes):
    Tree()
{
    if (!nodes.empty()){
        Point lower, upper;
        std::tie(lower, upper) = findVertices(nodes.cbegin(), nodes.cend());
        construct(_root, nodes.begin(), nodes.end(), lower, upper, 0);
        
        // The first throw is to check where construct() is implemented correctly. This throw right here might even be a bad design.
        //if (!hasValidBoxTypes(_root)) throw std::runtime_error("ERROR: At least one bounding box does not enclose its contents.");
        if (hasOverlappingContents(_root)) throw std::invalid_argument("ERROR: There are overlapping contents.");
        nodes.clear();
    }
}

void kdTree::insert(Node& node){
    if (_root->contentsOverlap(*node))
        throw std::invalid_argument("ERROR: The node to be inserted overlaps with the content of the tree.");
    // Find the deepest level that node can be inserted
}

Shape* kdTree::nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const{
    if (!_root){
        s = NAN;
        return nullptr;
    }
    if (current){
        s = current->distanceToSurface(pos, dir);
        if (current->encloses(pos) || (current->surfaceContains(pos) && s > 0.0)) return current;
    }

    s = std::numeric_limits<double>::max();
    std::stack<BoxType*> stack;
    stack.push(dynamic_cast<BoxType*>(_root.get()));
    Shape* next = nullptr;
    
    while (!stack.empty()){
        BoxType* box = stack.top();
        stack.pop();
        // Check if box has any viable candidates
        double dist = box->distanceToSurface(pos, dir);
        if (!box->encloses(pos) && !(box->surfaceContains(pos) && dist > 0.0)){
            // pos is outside of box, so the distance to any Shape inside of it will be at least dist
            if (s < dist) continue;
        }

        // Check the children, selectively mark those that the Point won't reach as visited
        for (unsigned i = 0; i < box->size(); i++){
            Shape* candidate = (*box)[i].get();
            dist = candidate->distanceToSurface(pos, dir);
            if (std::isnan(dist)) continue;

            if (BoxType* subbox = dynamic_cast<BoxType*>(candidate)) stack.push(subbox);
            else{
                if (candidate->encloses(pos)){
                    // Point is inside the Shape
                    s = dist;
                    return candidate;
                }

                if (candidate->surfaceContains(pos)){
                    if (dist == 0.0) continue; // Travelling away from candidate
                    if ((!current) || (current != candidate && current->surfaceContains(pos))){
                        // Current is nullptr or Point is on the intersection of current and candidate
                        if (!dir.isOrthogonal(candidate->normal(pos))){
                            // checks against tangency
                            s = 0.0;
                            return candidate;
                        }
                        continue;
                    }
                    if (current == candidate && dist > 0.0){
                        // Traveling within candidate
                        s = dist;
                        return candidate;
                    }
                }

                if (s > dist){
                    Point nextPos = pos;
                    nextPos.advance(dir, dist);
                    if (!dir.isOrthogonal(candidate->normal(nextPos))){
                        // checks against tangency
                        s = dist;
                        next = candidate;
                    }
                }
            }
        }
    }
    return next;
}

kdTree::iterator kdTree::medianNode(BoxType& bbox0, BoxType& bbox1, iterator begin, iterator end, const char axis){
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
    const auto N = end-begin;
    iterator median = begin+N/2;
    std::nth_element(begin, median, end, comp); // Find the median
    
    // bbox1 now should encloses the upper half of the nodes in nodes, now partialy sorted about element number N/2
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

void kdTree::construct(Node& current, iterator begin, iterator end, Point& lower, Point& upper, std::size_t level){
    /**
     * Inputs:
     *  current: pointer to the BoxType at which space partitioning occurs
     *  nodes: a vector of all nodes contained within *current
     *  begin, end: indices of nodes to partition
     *  lower, upper: vertices of the BoxType that must enclose all objects in nodes
     * Outputs:
     *  Recursively construct a kdTree as an octree
    **/

    if (!current){ // current is nullptr
        std::unique_ptr<BoxType> node = std::make_unique<BoxType>(lower, upper);
        node->setLevel(level);
        if (end - begin <= 8){
            for (auto it = begin; it != end; it++) (*node)[it-begin] = std::move(*it);
        } else{
            std::vector<BoxType> bbox(9);
            std::vector<iterator> sep(9, begin);
            bbox[0] = BoxType(lower, upper);
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

void kdTree::destruct(Node& current, NodeList& nodes){
    if (current){ // not nullptr
        if (BoxType* box = dynamic_cast<BoxType*>(current.get())){
            // If current is a BoxType, recursively destruct its children
            for (std::size_t i = 0; i < box->size(); i++) destruct((*box)[i], nodes);
        } else nodes.emplace_back(std::move(current));
        current.reset();
    }
}

bool kdTree::hasOverlappingContents(const Node& current) const{
    // Potentially O(N^2) in time complexity, 
    if (!current) return false; // nullptr;
    if (const BoxType* box = dynamic_cast<const BoxType*>(current.get())){
        for (std::size_t i = 0; i < box->size(); i++){
            if (!(*box)[i]) continue; // nullptr
            if (hasOverlappingContents((*box)[i])) return true;

            for (std::size_t j = i+1; j < box->size(); j++){
                if (!(*box)[j]) continue; // nullptr
                if ((*box)[i]->contentsOverlap(*(*box)[j])) return true;
            }
        }
    }
    return false; // not a BoxType
}

// bool hasValidBoxTypes(const Node& current){
//     if (!current) return true; // nullptr
//     if (const BoxType* box = dynamic_cast<const BoxType*>(current.get())){
//         for (std::size_t i = 0; i < box->size(); i++){
//             if (!(*box)[i]) continue; // nullptr
//             if (!box->encloses( *(*box)[i] )) return false;
//             if (!hasValidBoxTypes((*box)[i])) return false;
//         }
//         return true;
//     }
//     return true; // not a BoxType;
// }