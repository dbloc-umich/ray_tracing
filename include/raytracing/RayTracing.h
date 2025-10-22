#ifndef RAY_TRACING_H
#define RAY_TRACING_H

class Ray;
template<typename T> class Tree;

template<typename T>
double intensity (const Tree<T>& tree, const Ray& ray, bool isRefracted = true, bool isReflected = true) noexcept;

#endif