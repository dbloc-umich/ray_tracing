#ifndef RAY_TRACING_H
#define RAY_TRACING_H

class Direction;
class kDTree;
class Point;
class Shape;

Direction reflected(const Direction& in, const Direction& normal);
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2);

double intensity(const kDTree& tree, Point p, const Direction& dir, Shape* current=nullptr, double initial=1.0);

#endif // RAY_TRACING_H
