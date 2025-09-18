#include "Octree.h"
#include "Direction.h"
#include "DisjointBoundingBox.h"

#include <cmath>
#include <exception>
#include <stack>

namespace{
    using BoxType = DisjointBoundingBox;
}

Octree::Octree(NodeList& nodes):
    Tree()
{
    if (!nodes.empty()){
        Point lower, upper;
        std::tie(lower, upper) = findVertices(nodes.cbegin(), nodes.cend());
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

        // Check the contents first
        for (unsigned j = 0; j < box->numContents(); j++){
            Shape* candidate = (*box)(j).get(); // candidate does not point to a BoxType
            dist = candidate->distanceToSurface(pos, dir);
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

        // Check the children next, selectively mark those that the Point won't reach as visited
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


void Octree::construct(Node& current, NodeList& nodes, const Point& lower, const Point& upper, std::size_t level){
    if (!current){ // current is nullptr
        std::unique_ptr<BoxType> box = std::make_unique<BoxType>(lower, upper);
        box->setLevel(level);
        if (nodes.size() <= 8) {
            // There are 8 items or less, move them all into the children nodes
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

void Octree::destruct(Node& current, NodeList& nodes){
    if (current){ // not nullptr
        if (BoxType* box = dynamic_cast<BoxType*>(current.get())){
            // If current is a BoxType, recursively destruct its children
            for (std::size_t i = 0; i < box->size(); i++) destruct((*box)[i], nodes);
            for (std::size_t j = 0; j < box->numContents(); j++) nodes.emplace_back(std::move((*box)(j)));
        } else nodes.emplace_back(std::move(current));
        current.reset();
    }
}

bool Octree::hasOverlappingContents(const Node& current) const{
    // Potentially O(N^2) in time complexity

    if (!current) return false; // nullptr;
    if (const BoxType* box = dynamic_cast<const BoxType*>(current.get())){
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
    return false; // not a BoxType
}