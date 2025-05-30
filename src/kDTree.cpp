// #include "kDTree.h"
// #include "Direction.h"

// #include <algorithm>
// #include <cmath>
// #include <exception>
// #include <limits>
// #include <tuple>

// namespace{

// iterator medianNode(BoundingBox& bbox0, BoundingBox& bbox1, iterator begin, iterator end, const char axis='x'){
//     /**
//      * Inputs:
//      *  bbox0 is the bounding box that completely encloses all objects in nodes[begin:end)
//      *  bbox1 must be empty
//      * Outputs:
//      *  since this function is only called when N = end-begin > 8:
//      *      returns an iterator pointing to an element N/2 away from begin, where N=end-begin
//      *      the underlying nodeList is partially sorted such that the elements in [begin, median) are contained in bbox0, and those in [median, end) are contained in bbox1
//      *      bbox0 is now modified to only contain the lower half of nodes0
//     **/
    
//     auto comp = [axis](const Node& a, const Node& b)
//         {
//             if (axis == 'x') return a->xMin() < b->xMin();
//             else if (axis == 'y') return a->yMin() < b->yMin();
//             else return a->zMin() < b->zMin();
//         };
//     const auto N = end-begin;
//     iterator median = begin+N/2;
//     std::nth_element(begin, median, end, comp); // Find the median
    
//     // bbox1 now should encloses the upper half of the nodes in nodes, now partialy sorted about element number N/2
//     Point p = bbox0.lowerVertex();
//     if (axis == 'x') p.setX((*median)->xMin() - eps);
//     else if (axis == 'y') p.setY((*median)->yMin() - eps);
//     else p.setZ((*median)->zMin() - eps);
//     bbox1.setVertices(p, bbox0.upperVertex());

//     // Restructure bbox0 to only enclose the lower half of nodes0
//     Point lv, uv;
//     std::tie(lv, uv) = findVertices(begin, median);
//     bbox0.setVertices(lv, uv);
//     return median;
// }

// void construct(Node& current, iterator begin, iterator end, Point& lower, Point& upper, std::size_t level){
//     /**
//      * Inputs:
//      *  current: pointer to the BoundingBox at which space partitioning occurs
//      *  nodes: a vector of all nodes contained within *current
//      *  begin, end: indices of nodes to partition
//      *  lower, upper: vertices of the BoundingBox that must enclose all objects in nodes
//      * Outputs:
//      *  Recursively construct a kDTree as an octree
//     **/

//     if (!current){ // current is nullptr
//         std::unique_ptr<BoundingBox> node = std::make_unique<BoundingBox>(lower, upper);
//         node->setLevel(level);
//         if (end - begin <= 8){
//             for (auto it = begin; it != end; it++) (*node)[it-begin] = std::move(*it);
//         } else{
//             std::vector<BoundingBox> bbox(9);
//             std::vector<iterator> sep(9, begin);
//             bbox[0] = BoundingBox(lower, upper);
//             sep[8] = end;

//             // Divide a single box along the x-axis
//             sep[4] = medianNode(bbox[0], bbox[4], begin, end, 'x');

//             // Divide the 2 boxes along the y-axis
//             for (std::size_t i = 0; i < 8; i+=4){
//                 if (sep[i+4]-sep[i] > 8){
//                     sep[i+2] = medianNode(bbox[i], bbox[i+2], sep[i], sep[i+4], 'y');
//                 }
//             }

//             // Divide the 4 boxes along the z-axis
//             for (std::size_t i = 0; i < 8; i+=2){
//                 bool splitsInZ = (sep[i+2]-sep[i] > (int)8) && (sep[i+2]-sep[i] < end-begin);
//                 if (splitsInZ) sep[i+1] = medianNode(bbox[i], bbox[i+1], sep[i], sep[i+2], 'z');
//             }

//             for (std::size_t i = 1; i < 8; ){
//                 if (sep[i] == begin){
//                     bbox.erase(bbox.begin()+i);
//                     sep.erase(sep.begin()+i);
//                 } else i++;
//             }
//             bbox.pop_back();

//             for (std::size_t i = 0; i < bbox.size(); i++){
//                 Point lv = bbox[i].lowerVertex();
//                 Point uv = bbox[i].upperVertex();
//                 construct((*node)[i], sep[i], sep[i+1], lv, uv, level+1);
//             } 
//         }
//         current = std::move(node);
//     }    
// }

// void destruct(Node& current, NodeList& nodes){
//     if (current){ // not nullptr
//         if (BoundingBox* box = dynamic_cast<BoundingBox*>(current.get())){
//             // If current is a BoundingBox, recursively destruct its children
//             for (std::size_t i = 0; i < box->size(); i++) destruct((*box)[i], nodes);
//         } else nodes.emplace_back(std::move(current));
//         current.reset();
//     }
// }

// bool hasValidBoundingBoxes(const Node& current){
//     if (!current) return true; // nullptr
//     if (const BoundingBox* box = dynamic_cast<const BoundingBox*>(current.get())){
//         for (std::size_t i = 0; i < box->size(); i++){
//             if (!(*box)[i]) continue; // nullptr
//             if (!box->encloses( *(*box)[i] )) return false;
//             if (!hasValidBoundingBoxes((*box)[i])) return false;
//         }
//         return true;
//     }
//     return true; // not a BoundingBox;
// }

// bool hasOverlappingContents(const Node& current){
//     // Potentially O(N^2) in time complexity, 
//     if (!current) return false; // nullptr;
//     if (const BoundingBox* box = dynamic_cast<const BoundingBox*>(current.get())){
//         for (std::size_t i = 0; i < box->size(); i++){
//             if (!(*box)[i]) continue; // nullptr
//             if (hasOverlappingContents((*box)[i])) return true;

//             for (std::size_t j = i+1; j < box->size(); j++){
//                 if (!(*box)[j]) continue; // nullptr
//                 if ((*box)[i]->contentsOverlap(*(*box)[j])) return true;
//             }
//         }
//     }
//     return false; // not a BoundingBox
// }

// Shape* nextExternalNode(BoxStack& stack, const Point& pos, const Direction& dir, std::vector<Shape*>& visitedNodes, double& s){
//     if (stack.empty()) return nullptr;

//     BoundingBox* box = stack.top(); // stack.top() is always at least the parent of all nodes in visitedNodes
//     std::size_t index = box->max_size();
//     s = std::numeric_limits<double>::max();
//     for (std::size_t i = 0; i < box->size(); i++){
//         if (std::find(visitedNodes.cbegin(), visitedNodes.cend(), (*box)[i].get()) != visitedNodes.cend()) continue; // currently on a visited Node, will skip
//         double dist = (*box)[i]->distanceToSurface(pos, dir);
//         if (!std::isnan(dist)){
//             Point nextPos = pos;
//             nextPos.advance(dir, dist); // to check against tangency
//             if (!dir.isOrthogonal((*box)[i]->normal(nextPos)) && dist < s){
//                 index = i;
//                 s = dist;
//             }
//         }
//     }
    
//     if (index == box->max_size()){ // Can't find nearest neighbor at currentNode
//         if (visitedNodes.size() == stack.top()->size()){
//             visitedNodes.clear();
//             visitedNodes.push_back(stack.top());
//             stack.pop();
//         }
//         else if (visitedNodes.size() > 1) visitedNodes.push_back(box);
//         else{
//             if (!dynamic_cast<BoundingBox*>(visitedNodes[0])){
//                 stack.pop();
//                 visitedNodes[0] = box;
//             }
//             else visitedNodes.push_back(box);
//         }
//         return nextExternalNode(stack, pos, dir, visitedNodes, s);
//     }
//     if (BoundingBox* subbox = dynamic_cast<BoundingBox*>((*box)[index].get())){ // Goes deeper into currentNode
//         stack.push(subbox);
//         return nextExternalNode(stack, pos, dir, visitedNodes, s);
//     }
//     return (*box)[index].get(); // Found it
// }
// }

// kDTree::kDTree(NodeList& nodes):
//     Tree()
// {
//     if (!nodes.empty()){
//         Point lower, upper;
//         std::tie(lower, upper) = findVertices(nodes.cbegin(), nodes.cend());
//         construct(_root, nodes.begin(), nodes.end(), lower, upper, 0);
        
//         // The first throw is to check where construct() is implemented correctly. This throw right here might even be a bad design.
//         //if (!hasValidBoundingBoxes(_root)) throw std::runtime_error("ERROR: At least one bounding box does not enclose its contents.");
//         if (hasOverlappingContents(_root)) throw std::invalid_argument("ERROR: There are overlapping contents.");
//         nodes.clear();
//     }
// }

// void kDTree::insert(Node& node){
//     if (_root->contentsOverlap(*node))
//         throw std::invalid_argument("ERROR: The node to be inserted overlaps with the content of the tree.");
//     // Find the deepest level that node can be inserted
// }

// NodeList kDTree::remove(Node& node){
//     NodeList list;
//     destruct(node, list);
//     return list;
// }

// Shape* kDTree::nextNode(const Point& pos, const Direction& dir, Shape* current, double& s) const{
//     // Check current to see whether the Point is inside the Shape or it is entering from the surface
//     if (current){
//         s = current->distanceToSurface(pos, dir);
//         if (current->encloses(pos) || (current->surfaceContains(pos) && s > 0.0)) return current;
//     }

//     // current is guaranteed to not be the returned pointer
//     if (_root && !std::isnan(_root->distanceToSurface(pos, dir))){
//         BoxStack stack;
//         stack.emplace(dynamic_cast<BoundingBox*>(_root.get()));
//         std::vector<Shape*> visited{current};
//         visited.reserve(8);
//         auto nextNode = nextExternalNode(stack, pos, dir, visited, s);
//         if (!nextNode) s = NAN;
//         return nextNode;
//     }
//     s = NAN;
//     return nullptr;
// }