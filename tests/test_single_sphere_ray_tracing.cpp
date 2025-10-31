#include "catch2/catch2.hpp"

#include "ConstantProperty.h"
#include "Direction.h"
#include "kdTree.h"
#include "Material.h"
#include "Ray.h"
#include "RayTracing.h"
#include "Sphere.h"

#include <cmath>

TEST_CASE("one sphere, streamline with normal incidence"){
    double r = 2;
    double n = 1.3;
    double alpha = 0.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r, 0, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, HopMode::Streamline);
    double Itrue = exp(-2*alpha*r);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, refraction and reflection with normal incidence"){
    double r = 2;
    double n = 1.3;
    double alpha = 0.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r, 0, 0), Direction(1, 0, 0));
    double A = exp(-2*alpha*r);
    double R = (1.0-n)/(1.0+n); R *= R;
    double T = 1.0 - R;
    double I = intensity(tree, ray, HopMode::Detailed);
    double Itrue = R + T*T*A/(1-R*A);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, streamline with non-normal incidence"){
    double r = 2;
    double n = 1.3;
    double alpha = 0.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r, 0.5*r, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, HopMode::Streamline);
    double Itrue = exp(-std::sqrt(3)*alpha*r);
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, refraction and reflection with non-normal incidence"){
    double r = 2;
    double n = 1.3;
    double alpha = 0.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, mat));
    kdTree tree(nodes);

    Ray ray(Point(-2*r, 0.5*r, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, HopMode::Detailed);
    double Itrue = 0.338424;
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}

TEST_CASE("one sphere, tangential incidence"){
    double r = 2;
    double n = 1.3;
    double alpha = 0.3;
    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(n));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(alpha));
    std::vector<std::unique_ptr<Shape>> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, mat));
    kdTree tree(nodes);

    Ray ray(Point(-r, r, 0), Direction(1, 0, 0));
    double I = intensity(tree, ray, HopMode::Detailed);
    double Itrue = 1.0;
    REQUIRE(fabs(I-Itrue)/Itrue < 1e-6);
}