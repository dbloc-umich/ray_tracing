#include "kdTree.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <stack>

#include "Box.h"
#include "Direction.h"
#include "Point.h"

kdTree::kdTree(PtrList& ptrs):
    Tree()
{
    if (!ptrs.empty()){
        _root = Node(boundingBox(ptrs.cbegin(), ptrs.cend()));
        construct(_root, ptrs.begin(), ptrs.end(), 0, 'x');
        
        ptrs.clear();
        if (hasOverlappingContents(_root)){
            // Returns the shapes back to the original container, but the correct order is not guaranteed;
            destruct(_root, ptrs);
            throw std::invalid_argument("ERROR: There are overlapping contents.");
        }  
    }
}

void kdTree::insert(std::unique_ptr<Shape> shape){
    // Empty tree, make a box over the shape
    if (!_root){
        PtrList ptrs;
        ptrs.push_back(std::move(shape));
        _root = Node(boundingBox(ptrs.cbegin(), ptrs.cend()));
        _root[0] = std::make_unique<Node>(std::move(ptrs[0]));
    }
    
    // Collision check
    if (_root.leavesOverlap(*shape)){
        throw std::invalid_argument("ERROR: The node to be inserted overlaps with the content of the tree.");
    }
    
    PtrList ptrs;
    ptrs.push_back(std::move(shape));
    // Check if the Shape can be put inside _root
    if (_root->encloses(*ptrs[0])){
        Node& node = smallestBox(_root, ptrs[0]); // the deepest level that node can be inserted
        // Check if there's still room to insert without destroying the subtree
        for (std::size_t i = 0; i < 2; i++){
            if (!node[i]){
                node[i] = std::make_unique<Node>(std::move(ptrs[0]));
                return;
            }
        }
        // Node is full, now building a brand new subtree
        std::size_t level = node.level();
        char axis = level%3 == 0 ? 'x' : (level%3 == 1 ? 'y' : 'z');
        destruct(node, ptrs);
        node = Node(boundingBox(ptrs.cbegin(), ptrs.cend()));
        construct(node, ptrs.begin(), ptrs.end(), level, axis);
    } else{
        destruct(_root, ptrs);
        _root = Node(boundingBox(ptrs.cbegin(), ptrs.cend()));
        construct(_root, ptrs.begin(), ptrs.end(), 0, 'x');
    }
}

Shape* kdTree::nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const noexcept{
    if (!_root){
        s = NAN;
        return nullptr;
    }
    if (current){
        s = current->distanceToSurface(pos, dir);
        if (current->encloses(pos) || (current->surfaceContains(pos) && s > 0.0)) return current;
    }

    s = std::numeric_limits<double>::max();
    std::stack<std::reference_wrapper<const Node>> stack; // keeps tracks of nodes to visit
    stack.emplace(_root);
    Shape* next = nullptr;

    while (!stack.empty()){
        auto& node = stack.top().get();
        stack.pop();
        // Check if box has any viable candidates
        double dist = node->distanceToSurface(pos, dir);
        if (!node->encloses(pos) && !(node->surfaceContains(pos) && dist > 0.0)){
            // pos is outside of this Shape, so the distance to any Shape inside of it will be at least dist
            if (s < dist) continue;
        }

        // Check the children, selectively mark those that the Point won't reach as visited
        for (std::size_t i = 0; i < node.max_size(); i++){
            if (!node[i]) continue;
            auto& candidate = *(node[i]);
            dist = candidate->distanceToSurface(pos, dir);
            if (std::isnan(dist)) continue;

            if (!candidate.isLeaf()) stack.emplace(candidate); // Not a leaf node, keeps adding to the stack
            else{
                // Reaches a leaf node, now performing distance calcs
                if (candidate->encloses(pos)){
                    // Point is inside the Shape
                    s = dist;
                    return candidate.get();
                }

                if (candidate->surfaceContains(pos)){
                    if (dist == 0.0) continue; // Travelling away from candidate
                    if ((!current) || (current != candidate.get() && current->surfaceContains(pos))){
                        // Current is nullptr or Point is on the intersection of current and candidate
                        if (!dir.isOrthogonal(candidate->normal(pos))){
                            // checks against tangency
                            s = 0.0;
                            return candidate.get();
                        }
                        continue;
                    }
                    if (current == candidate.get() && dist > 0.0){
                        // Traveling within candidate
                        s = dist;
                        return current;
                    }
                }

                if (s > dist){
                    Point nextPos = pos;
                    nextPos.advance(dir, dist);
                    if (!dir.isOrthogonal(candidate->normal(nextPos))){
                        // checks against tangency
                        s = dist;
                        next = candidate.get();
                    }
                }
            }
        }
    }
    return next;
}

kdTree::iterator kdTree::medianNode(Box& box0, Box& box1, iterator begin, iterator end, char axis) noexcept{
    /**
     * Inputs:
     *  box0 is the bounding box that completely encloses all objects in nodes[begin:end)
     *  box1 is a placeholder box
     * Outputs:
     *  since this function is only called when N = end-begin > 8:
     *      returns an iterator pointing to an element N/2 away from begin, where N=end-begin
     *      the underlying PtrList is partially sorted such that the elements in [begin, median) are contained in box0, and those in [median, end) are contained in box1
     *      box0 is now modified to only contain the lower half of nodes0
    **/
    
    auto comp = [axis](const std::unique_ptr<Shape>& a, const std::unique_ptr<Shape>& b)
        {
            if (axis == 'x') return a->xMin() < b->xMin();
            else if (axis == 'y') return a->yMin() < b->yMin();
            else return a->zMin() < b->zMin();
        };
    const auto N = end-begin;
    iterator median = begin+N/2;
    std::nth_element(begin, median, end, comp); // Find the median
    
    // box1 now should encloses the upper half of the nodes, now partialy sorted about element number N/2
    Point p = box0.lowerVertex();
    if (axis == 'x') p.setX((*median)->xMin() - eps);
    else if (axis == 'y') p.setY((*median)->yMin() - eps);
    else p.setZ((*median)->zMin() - eps);
    box1.setVertices(p, box0.upperVertex());

    // Restructure box0 to only enclose the lower half of nodes0
    Point uv = box0.upperVertex();
    if (axis == 'x') uv.setX(box0.xMin());
    else if (axis == 'y') uv.setY(box0.yMin());
    else uv.setZ(box0.zMin());

    for (auto& it = begin; it != median; it++){
        if (axis == 'x' && (*it)->xMax() > uv.x()) uv.setX((*it)->xMax());
        else if (axis == 'y' && (*it)->yMax() > uv.y()) uv.setY((*it)->yMax());
        else if (axis == 'z' && (*it)->zMax() > uv.z()) uv.setZ((*it)->zMax());
    }
    box0.setUpperVertex(uv);
    return median;
}

void kdTree::construct(Node& current, iterator begin, iterator end, std::size_t level, char axis) noexcept{
    /**
     * Inputs:
     *  current: Node pointing to the Box at which space partitioning occurs
     *  begin, end: pointers to the beginning and end of a PtrList containing all ptrs that fit inside box
     *  box: a Box object that must enclose all objects in [begin, end)
     * Outputs:
     *  Recursively construct a kdTree
    **/

    if (end - begin <= 2){ // 1 or 2 pointers left
        for (auto it = begin; it != end; it++){
            current[it-begin] = std::make_unique<Node>(std::move(*it));
        }
    } else{
        Box box0 = dynamic_cast<Box&>(*current);
        Box box1;
        iterator median = medianNode(box0, box1, begin, end, axis);
        current[0] = std::make_unique<Node>(std::make_unique<Box>(box0));
        current[1] = std::make_unique<Node>(std::make_unique<Box>(box1));

        char next;
        if (axis == 'x') next = 'y';
        else if (axis == 'y') next = 'z';
        else next = 'x';
        construct(*current[0], begin, median, level+1, next);
        construct(*current[1], median, end, level+1, next);
    }
}

void kdTree::destruct(Node& current, PtrList& ptrs) noexcept{
    // current must be non-null
    if (current.isLeaf()) ptrs.emplace_back(current.release());
    else{
        if (current[0]) destruct(*current[0], ptrs);
        if (current[1]) destruct(*current[1], ptrs);
        current.reset();
    }
}

Node& kdTree::smallestBox(Node& current, const std::unique_ptr<Shape>& shape) const noexcept{
    // current must be non-null and current->encloses(*shape) must be true
    for (std::size_t i = 0; i < 2; i++){
        if (current[i] && !current[i]->isLeaf() && (*current[i])->encloses(*shape)){
            return smallestBox(*current[i], shape);
        }
    }
    return current;
}

bool kdTree::hasOverlappingContents(const Node& current) const noexcept{
    // Should be O(N*logN) in time complexity
    if (!current || current.isLeaf()) return false; // nullptr;

    if (current[0] && current[1]){ // no null children
        return hasOverlappingContents(*current[0])
            || hasOverlappingContents(*current[1])
            || current[0]->leavesOverlap(*current[1]);
    }
    if (current[0]) return hasOverlappingContents(*current[0]);
    return hasOverlappingContents(*current[1]);
}