#include <Tree.h>
#include <cmath>

template<typename T>
double Tree<T>::xMin() const noexcept{ return _root ? _root->xMin() : NAN; }

template<typename T>
double Tree<T>::xMax() const noexcept{ return _root ? _root->xMax() : NAN; }

template<typename T>
double Tree<T>::yMin() const noexcept{ return _root ? _root->yMin() : NAN; }

template<typename T>
double Tree<T>::yMax() const noexcept{ return _root ? _root->yMax() : NAN; }

template<typename T>
double Tree<T>::zMin() const noexcept{ return _root ? _root->zMin() : NAN; }

template<typename T>
double Tree<T>::zMax() const noexcept{ return _root ? _root->zMax() : NAN; }

template<typename T>
std::vector<Node> Tree<T>::remove(Node& node){
    NodeList list;
    destruct(node, list);
    return list;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Tree<T>& tree){
    if (tree._root) os << *(tree._root);
    else os << "Empty tree.";
    return os;
}

template<typename T>
std::pair<Point, Point> Tree<T>::findVertices(const_iterator begin, const_iterator end) const{
    double x0, y0, z0, x1, y1, z1;
    x0 = y0 = z0 = std::numeric_limits<double>::max();
    x1 = y1 = z1 = std::numeric_limits<double>::min();
    for (const_iterator it = begin; it != end; it++){
        x0 = std::min(x0, (*it)->xMin());
        y0 = std::min(y0, (*it)->yMin());
        z0 = std::min(z0, (*it)->zMin());
        x1 = std::max(x1, (*it)->xMax());
        y1 = std::max(y1, (*it)->yMax());
        z1 = std::max(z1, (*it)->zMax());
    }
    const double epsx = eps*(x1-x0);
    const double epsy = eps*(y1-y0);
    const double epsz = eps*(z1-y0);
    x0 -= epsx; y0 -= epsy; z0 -= epsz;
    x1 += epsx; y1 += epsy; z1 += epsz;
    return std::make_pair(Point(x0, y0, z0), Point(x1, y1, z1));
}