#ifndef RAY_TRACING_H
#define RAY_TRACING_H

class Direction;
class Point;
class Shape;
class Tree;

Direction reflected(const Direction& in, const Direction& normal);
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2);

double intensity(const Tree& tree, Point p, const Direction& dir,
                 bool isRefracted = true, bool isReflected = true,
                 Shape* current = nullptr, double initial = 1.0);

#endif // RAY_TRACING_H
  