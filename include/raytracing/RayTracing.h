#ifndef RAY_TRACING_H
#define RAY_TRACING_H

class Direction;
class Octree;
class Point;
class Shape;

Direction reflected(const Direction& in, const Direction& normal);
Direction refracted(const Direction& in, const Direction& normal, double n1, double n2);

double intensity(const Octree& tree, Point p, const Direction& dir,
                 bool isRefracted = false, bool isReflected = false,
                 Shape* current = nullptr, double initial = 1.0);

#endif // RAY_TRACING_H
