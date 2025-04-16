#include "Tree.h"
#include <cmath>

decltype(auto) findVertices(Tree::const_iterator cbegin, Tree::const_iterator cend){
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


double Tree::xMin() const noexcept{ return _root ? _root->xMin() : NAN; }
double Tree::xMax() const noexcept{ return _root ? _root->xMax() : NAN; }
double Tree::yMin() const noexcept{ return _root ? _root->yMin() : NAN; }
double Tree::yMax() const noexcept{ return _root ? _root->yMax() : NAN; }
double Tree::zMin() const noexcept{ return _root ? _root->zMin() : NAN; }
double Tree::zMax() const noexcept{ return _root ? _root->zMax() : NAN; }

std::ostream& operator<<(std::ostream& os, const Tree& tree){
    if (tree._root) os << *(tree._root);
    else os << "Empty tree.";
    return os;
}