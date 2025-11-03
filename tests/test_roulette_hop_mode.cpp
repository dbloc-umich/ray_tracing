// Testing whether the aggregate result from HopMode::Roulette differs significantly from HopMode::Detailed

#include "catch2/catch2.hpp"

#include "ConstantProperty.h"
#include "Direction.h"
#include "kdTree.h"
#include "Material.h"
#include "Ray.h"
#include "RayTracing.h"
#include "Sphere.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <random>

TEST_CASE("roulette_test"){
    double L = 10;
    double V0 = 8*L*L*L;
    double pf = 0.1;
    std::size_t N = 100; // number of droplets
    std::size_t Ny = 200; // number of rays on one axis

    std::shared_ptr<Material> mat = std::make_shared<Material>();
    mat->addProperty(Prop::refractiveIndex, std::make_unique<ConstantProperty>(1.3));
    mat->addProperty(Prop::attenuationCoefficient, std::make_unique<ConstantProperty>(0.3));

    // Create a random distribution of trees
    double r = std::cbrt(V0*pf/N/(4.0*M_PI/3));
    kdTree tree;
    std::default_random_engine rng(1);
    std::uniform_real_distribution<double> xDist(-(L-r), L-r);
    std::uniform_real_distribution<double> yDist(-(L-r), L-r);
    std::uniform_real_distribution<double> zDist(-(L-r), L-r);

    unsigned i = 0;
    while (i < N){
        Point r0(xDist(rng), yDist(rng), zDist(rng));
        std::unique_ptr<Shape> sph = std::make_unique<Sphere>(r0, r, mat);
        bool overlaps = tree.leavesOverlap(*sph); // Check for collision
        if (!overlaps){
            tree.insert(std::move(sph));
            i++;
        }
    }

     // Ray tracing
    const Direction dir(1, 0, 0);
    double I1 = 0.0, I2 = 0.0;
    double x = tree.xMin() - L;
    double dy = (tree.yMax() - tree.yMin())/Ny;
    double dz = (tree.zMax() - tree.zMin())/Ny;

    for (unsigned iy = 0; iy < Ny+1; iy++){
        double y = tree.root()->yMin() + dy*iy;
        for (unsigned iz = 0; iz < Ny+1; iz++){
            double z = tree.root()->zMin() + dz*iz;                       
            Ray ray(Point(x, y, z), dir);
            I1 += intensity(tree, ray, HopMode::Roulette);

            ray = Ray(Point(x, y, z), dir);
            I2 += intensity(tree, ray, HopMode::Detailed);
        }
    }
    REQUIRE(fabs(I1-I2)/I2 < 1e-3);
}