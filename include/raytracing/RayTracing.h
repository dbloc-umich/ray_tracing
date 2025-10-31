#ifndef RAY_TRACING_H
#define RAY_TRACING_H

class Ray;
template<typename T> class Tree;

/**
 * The enum class HopMode is used to indicate how a ray should be treated upon entering an interface.
 * There are three possible values:
 *  Streamline: The ray penetrates through the interface without changing direction or losing intensity.
 *              This is a simpliest treatment, but the most inaccurate.
 *  Roulette: The ray stochastically choses if it's reflected or transmitted, without losing intensity.
 *  Detailed: The ray is split into a reflected ray and possibly a transmitted ray, each with its own direction
 *            and intensity. Not recommended for large problems.
**/
enum class HopMode{ Streamline, Roulette, Detailed };

template<typename T>
double intensity (const Tree<T>& tree, Ray& ray, const HopMode& mode = HopMode::Roulette) noexcept;

#endif