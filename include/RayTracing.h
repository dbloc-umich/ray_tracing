#ifndef RAY_TRACING_H
#define RAY_TRACING_H

class Direction;
class kDTree;
class Point;
class Shape;

Direction reflected(const Direction& in, const Direction& normal);
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2);

double intensity(const kDTree& tree, Point p, const Direction& dir, Shape* current=nullptr, double initial=1.0);
//std::pair<double, double> p_leak(const Octree& tree, const Direction& dir);

#endif // RAY_TRACING_H
