#include "catch2/catch2.hpp"
#include "Ray.h"
#include "Direction.h"
#include "Sphere.h"
#include "Octree.h"

#include <cmath>

TEST_CASE("one sphere, streamline with normal incidence"){
    double r = 2;
    double Sigma_t = 0.3;
    double n = 1.3;
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, n));
    Octree tree(nodes);

    Ray ray(Point(-2*r, 0, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, false, false);
    double Itrue = exp(-2*Sigma_t*r);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, refraction and reflection with normal incidence"){
    double r = 2;
    double Sigma_t = 0.3;
    double n = 1.3;
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, n));
    Octree tree(nodes);

    Ray ray(Point(-2*r, 0, 0), Direction(1, 0, 0));
    double A = exp(-2*Sigma_t*r);
    double R = (1.0-n)/(1.0+n); R *= R;
    double T = 1.0 - R;
    double I = intensity(tree, ray, true, true);
    double Itrue = R + T*T*A/(1-R*A);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, streamline with non-normal incidence"){
    double r = 2;
    double Sigma_t = 0.3;
    double n = 1.3;
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, n));
    Octree tree(nodes);

    Ray ray(Point(-2*r, 0.5*r, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, false, false);
    double Itrue = exp(-std::sqrt(3)*Sigma_t*r);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, refraction and reflection with non-normal incidence"){
    double r = 2;
    double Sigma_t = 0.3;
    double n = 1.3;
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, n));
    Octree tree(nodes);

    Ray ray(Point(-2*r, 0.5*r, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, true, true);
    double Itrue = 0.338424;
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, tangential incidence"){
    double r = 2;
    double Sigma_t = 0.3;
    double n = 1.3;
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, n));
    Octree tree(nodes);

    Ray ray(Point(-r, r, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, true, true);
    double Itrue = 1.0;
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}