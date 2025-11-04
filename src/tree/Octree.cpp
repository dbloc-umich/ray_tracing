// #include "Octree.h"

// #include <cmath>
// #include <limits>
// #include <stack>

// #include "Box.h"
// #include "Direction.h"
// #include "Point.h"
// #include "Shape.h"

// Octree::Octree(PtrList& ptrs):
//     Tree()
// {
//     if (!ptrs.empty()){
//         _root = UNode(boundingBox(ptrs.cbegin(), ptrs.cend()));
//         construct(_root, ptrs);
        
//         if (hasOverlappingContents(_root)){
//             // Returns the shapes back to the original container, but the correct order is not guaranteed;
//             destruct(_root, ptrs);
//             throw std::invalid_argument("ERROR: There are overlapping contents.");
//         }  
//     }
// }

// void Octree::insert(std::unique_ptr<Shape> shape){
//     // Empty tree, make a box over the shape
//     if (!_root){
//         PtrList ptrs;
//         ptrs.push_back(std::move(shape));
//         _root = UNode(boundingBox(ptrs.cbegin(), ptrs.cend()));
//         _root[0] = std::make_unique<UNode>(std::move(ptrs[0]));
//     } else{
//         // Collision check
//         if (_root.leavesOverlap(*shape)){
//             throw std::invalid_argument("ERROR: The node to be inserted overlaps with the content of the tree.");
//         }
        
//         PtrList ptrs;
//         ptrs.push_back(std::move(shape));
//         // Check if the Shape can be put inside _root
//         if (_root->encloses(*ptrs[0])){
//             UNode& node = smallestBox(_root, ptrs[0]); // the deepest level that node can be inserted

//             if (node.size() != node.max_size()){
//                 // Node is not yet full, may be able to put the shape directly under _children
//                 // Check if all the existing children are leaf nodes
//                 bool isLastInternal = true;
//                 for (std::size_t oct = 0; oct < node.max_size(); oct++){
//                     if (node[oct] && !node[oct]->isLeaf()){
//                         isLastInternal = false;
//                         break;
//                     }
//                 }
//                 if (isLastInternal){
//                     // insert directly to the first empty children
//                     for (std::size_t oct = 0; oct < node.max_size(); oct++){
//                         if (!node[oct]){
//                             node[oct] = std::make_unique<UNode>(std::move(ptrs[0]));
//                             break;
//                         }
//                     }
//                 }
//             }
            
//             // Check if the Shape can be added to the node's direct leaves
//             std::size_t oct = node.octant(*ptrs[0]);
//             if (oct == 8) node.emplace(std::make_unique<Node>(std::move(ptrs[0])));
//             else{
//                 if (!node[oct]) node[oct] = std::make_unique<UNode>(std::move(ptrs[0]));
//                 else{
//                     // destruct the children but leave the leaves alone, then reconstruct the current node
//                     for (oct = 0; oct < node.max_size(); oct++){
//                         if (node[oct]) destruct(*node[oct], ptrs);
//                     }
//                     construct(node, ptrs, node.level());
//                 }
//             }
//         } else{
//             destruct(_root, ptrs);
//             _root = UNode(boundingBox(ptrs.cbegin(), ptrs.cend()));
//             construct(_root, ptrs, 0);
//         }
//     }
// }

// Shape* Octree::nextShape(const Point& pos, const Direction& dir, Shape* current, double& s) const noexcept{
//     if (!_root){
//         s = NAN;
//         return nullptr;
//     }
//     if (current){
//         s = current->distanceToSurface(pos, dir);
//         if (current->encloses(pos) || (current->surfaceContains(pos) && s > 0.0)) return current;
//     }

//     s = std::numeric_limits<double>::max();
//     std::stack<std::reference_wrapper<const UNode>> stack; // keeps tracks of nodes to visit
//     stack.emplace(_root);
//     Shape* next = nullptr;
    
//     while (!stack.empty()){
//         auto& node = stack.top().get();
//         stack.pop();
//         // Check if box has any viable candidates
//         double dist = node->distanceToSurface(pos, dir);
//         if (!node->encloses(pos) && !(node->surfaceContains(pos) && dist > 0.0)){
//             // pos is outside of this Shape, so the distance to any Shape inside of it will be at least dist
//             if (s < dist) continue;
//         }

//         // Check the leaves
//         for (std::size_t j = 0; j < node.leafCount(); j++){
//             auto& candidate = *(node(j));
//             dist = candidate->distanceToSurface(pos, dir);
//             if (std::isnan(dist)) continue;
//             if (candidate->encloses(pos)){
//                 // Point is inside the Shape
//                 s = dist;
//                 return candidate.get();
//             }

//             if (candidate->surfaceContains(pos)){
//                 if (dist == 0.0) continue; // Travelling away from candidate
//                 if ((!current) || (current != candidate.get() && current->surfaceContains(pos))){
//                     // Current is nullptr or Point is on the intersection of current and candidate
//                     if (!dir.isOrthogonal(candidate->normal(pos))){
//                         // checks against tangency
//                         s = 0.0;
//                         return candidate.get();
//                     }
//                     continue;
//                 }
//                 if (current == candidate.get() && dist > 0.0){
//                     // Traveling within candidate
//                     s = dist;
//                     return current;
//                 }
//             }

//             if (s > dist){
//                 Point nextPos = pos;
//                 nextPos.advance(dir, dist);
//                 if (!dir.isOrthogonal(candidate->normal(nextPos))){
//                     // checks against tangency
//                     s = dist;
//                     next = candidate.get();
//                 }
//             }
//         }

//         // Check the children, selectively mark those that the Point won't reach as visited
//         for (std::size_t i = 0; i < node.max_size(); i++){
//             if (!node[i]) continue;
//             auto& candidate = *(node[i]);
//             dist = candidate->distanceToSurface(pos, dir);
//             if (std::isnan(dist)) continue;

//             if (!candidate.isLeaf()) stack.emplace(candidate); // Not a leaf node, keeps adding to the stack
//             else{
//                 // Reaches a leaf node, now performing distance calcs
//                 if (candidate->encloses(pos)){
//                     // Point is inside the Shape
//                     s = dist;
//                     return candidate.get();
//                 }

//                 if (candidate->surfaceContains(pos)){
//                     if (dist == 0.0) continue; // Travelling away from candidate
//                     if ((!current) || (current != candidate.get() && current->surfaceContains(pos))){
//                         // Current is nullptr or Point is on the intersection of current and candidate
//                         if (!dir.isOrthogonal(candidate->normal(pos))){
//                             // checks against tangency
//                             s = 0.0;
//                             return candidate.get();
//                         }
//                         continue;
//                     }
//                     if (current == candidate.get() && dist > 0.0){
//                         // Traveling within candidate
//                         s = dist;
//                         return current;
//                     }
//                 }

//                 if (s > dist){
//                     Point nextPos = pos;
//                     nextPos.advance(dir, dist);
//                     if (!dir.isOrthogonal(candidate->normal(nextPos))){
//                         // checks against tangency
//                         s = dist;
//                         next = candidate.get();
//                     }
//                 }
//             }
//         }
//     }
//     return next;
// }


// void Octree::construct(UNode& current, PtrList& ptrs, std::size_t level) noexcept{
//     current.setLevel(level);
//     if (ptrs.size() <= 8) {
//         // There are 8 items or less, move them all into the children ptrs
//         for (std::size_t i = 0; i < ptrs.size(); i++){
//             current[i] = std::make_unique<UNode>(std::move(ptrs[i]));
//         }
//     }
//     else{
//         std::array<PtrList, 8> octantList;
//         while (!ptrs.empty()){
//             const std::size_t oct = current.octant(*(ptrs.back()));
//             if (oct == 8) current.emplace(std::make_unique<Node>(ptrs.back().release()));
//             else octantList[oct].emplace_back(ptrs.back().release());
//             ptrs.pop_back();
//         }

//         for (std::size_t oct = 0; oct < 8; oct++){
//             if (!octantList[oct].empty()){
//                 Point subLower(current->xMin(), current->yMin(), current->zMin());
//                 Point subUpper(current->xMax(), current->yMax(), current->zMax());
//                 double xMid = (current->xMin() + current->xMax())/2;
//                 double yMid = (current->yMin() + current->yMax())/2;
//                 double zMid = (current->zMin() + current->zMax())/2;

//                 if (oct < 4) subUpper.setX(xMid);
//                 else subLower.setX(xMid);
//                 if (oct % 4 < 2) subUpper.setY(yMid);
//                 else subLower.setY(yMid);
//                 if (!(oct % 2)) subUpper.setZ(zMid);
//                 else subLower.setZ(zMid);
                
//                 current[oct] = std::make_unique<UNode>(std::make_unique<Box>(subLower, subUpper));
//                 construct(*current[oct], octantList[oct], level+1);
//             }
//         }
//     }
// }    

// void Octree::destruct(UNode& current, PtrList& ptrs) noexcept{
//     // current must be non-null
//     if (current.isLeaf()) ptrs.emplace_back(current.release());
//     else{
//         for (std::size_t j = 0; j < current.leafCount(); j++){
//             ptrs.emplace_back(current(j)->release());
//             current(j).reset();
//         }
//         for (std::size_t i = 0; i < current.max_size(); i++){
//             if (current[i]) destruct(*current[i], ptrs);
//         }
//         current.reset();
//     }
// }

// UNode& Octree::smallestBox(UNode& current, const std::unique_ptr<Shape>& shape) const noexcept{
//     // current must be non-null and current->encloses(*shape) must be true
//     std::size_t oct = current.octant(*shape);
//     if (oct == 8 || !current[oct] || current[oct]->isLeaf()) return current;
//     return smallestBox(*current[oct], shape);
// }

// bool Octree::hasOverlappingContents(const UNode& current) const noexcept{
//     // Should be O(N*logN) at best, but the leaves can worsen it
//     if (!current || current.isLeaf()) return false; // nullptr;

//     for (std::size_t j = 0; j < current.leafCount(); j++){
//         for (std::size_t k = j+1; k < current.leafCount(); k++){
//             if ((*current(j))->overlaps(**current(k))) return true;
//         }
//     }

//     for (std::size_t i = 0; i < current.max_size(); i++){
//         if (!current[i]) continue;
//         if (hasOverlappingContents(*current[i])) return true;
//         for (std::size_t j = 0; j < current.max_size(); j++){
//             if (!current[j]) continue;
//             if (current[i]->leavesOverlap(*current[j])) return true;
//         }
//     }

//     return false;
// }