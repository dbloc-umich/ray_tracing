#include "catch2/catch2.hpp"

#include "ConstantProperty.h"
#include "Direction.h"
#include "kdTree.h"
#include "Material.h"
#include "Ray.h"
#include "RayTracing.h"
#include "Sphere.h"

#include <cmath>

TEST_CASE("two spheres, streamline with normal incidence"){
    double r1 = 2.5;
    double r2 = 1.5;
    double alpha = 0.3;
    double n = 1.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, mat));
    nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r1, 0, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, false, false);
    double Itrue = exp(-2*alpha*(r1+r2));
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("two spheres, refraction and reflection with normal incidence"){
    double r1 = 2.5;
    double r2 = 1.5;
    double alpha = 0.3;
    double n = 1.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, mat));
    nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r1, 0, 0), Direction(1, 0, 0));
    double A = exp(-2*alpha*(r1+r2));
    double R = (1-n)/(1+n); R *= R;
    double T = 1-R;
    double I = intensity(tree, ray, true, true);
    double Itrue = R + T*T*A/(1-R*A);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("two spheres, refraction and reflection with non-normal incidence"){
    double r1 = 2.5;
    double r2 = 1.5;
    double alpha = 0.3;
    double n = 1.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, mat));
    nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r1, 0.5*r1, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, true, true);
    double Itrue = 0.115967;
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("two spheres, tangential incidence to both spheres"){
    double r1 = 2.5;
    double r2 = 1.5;
    double alpha = 0.3;
    double n = 1.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, mat));
    nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, mat));
    kdTree tree(nodes);

    Ray ray(Point(r1, 0.5*r1, 0), Direction(0, -1, 0));
    double I = intensity(tree, ray, true, true);
    double Itrue = 1.0;
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}